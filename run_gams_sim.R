
# This script runs the GAM analysis on a SNP-by-SNP basis. Will loop through a set number of SNPs and collate results to make it possible to batch the analysis on a computing cluster.

# Read in command line arguments
	args<-commandArgs(trailingOnly=TRUE)

  print(args)

	task_id <- as.numeric(args[1])
  interval <- as.numeric(args[2])
  climate_var_dif <- as.character(args[3])
  data_file <- as.character(args[4])

# Load libraries
  library(mgcv, lib.loc = "/u/home/l/lukembro/R/x86_64-pc-linux-gnu-library/3.5") # Make sure to load latest version of mgcv
  library(tidyverse)

# Load data
  load(paste0("../", data_file))
 
## Read in snps to skip
 # skip_these_SNPs <- read_csv("../skip_these_SNPs_2018_12_22.csv")

# Save functions
# function to transform from raw to scaled variables
  forward_transform <- function(x, var, means, sds){
    ( x - means[var]) / sds[var]
  }
  
  # Function to back transform variables
  back_transform <- function(x, var, means, sds){
    (x * sds[var]) + means[var]
  }

 # Initialize variables
  gam_list <- list()
  x = 1
  
 # Convert fixed and random parts of formula to actual formula objetc
    make_formula <- function(fixed, random){
      return(formula(formula(enquote(paste0(fixed, random)))))
    }
    
  # Formula for random effects
    random_effects <-  '+ s(accession, bs = "re")'
  
  # Loop through snps by index based on task id
  for(snp_index in task_id:(task_id + interval - 1)){
        
    # To avoid going past number of snps
    if(snp_index > length(snp_col_names)){
      next
    }
    
    # For timing loops
    start_time <- Sys.time()
    
    # Choose snp
    snp <- snp_col_names[snp_index]
    
    # Skip SNPs that are already run
    # if(snp %in% skip_these_SNPs$snp){
    #   cat("Skipping SNP:", snp, " ... \n")
    #   next
    # }
    
   # Make sure there's at least 2 levels
    if(length(table(dat_snp_all_unscaled[, snp])) < 2){
      cat("Skipping because not at least 2 levels")
      next
    }

# Run gams
   	cat("Working on: ", snp, "... number: ", x, " ...\n" )

   ## Testing out residual gams
   	fixed_effects_resids <- paste0(paste0("rgr_resids ~ ", snp, " + s(tmax_dif, by = ", snp, ", bs=\"cr\")"))
   	
  # Gam with interaction effect
    gam_snp_int = bam(
                   formula = formula(fixed_effects_resids),
                  data = dat_snp_all_unscaled,
                  discrete = FALSE,
                  nthreads = 8,
                  method = "fREML")
                
    print(summary(gam_snp_int))
  
  # Save deviance differences and Rsquared into model
    gam_snp_int_summary <- summary(gam_snp_int)

    gam_snp_int$dev_dif <- NA
    gam_snp_int$rsq_dif <- NA
    
   # Save sample size per genotype
     gen_table <-  table(dat_snp_all_unscaled[, snp])
     gam_snp_int$n_gen0 <- as.numeric(gen_table["0"])
     gam_snp_int$n_gen1 <- as.numeric(gen_table["1"])
     gam_snp_int$n_gen2 <- as.numeric(gen_table["2"])
    

  # Make predictions

   # When SNPs are continuous, make own genotypes - Dont scale geontype
      gen_temp <- data.frame(accession = "1", 
                 genotype = c(0, 1, 2))            

      pred_sub <- left_join(pred, gen_temp)
   
  ## For SNPs as continuous
      pred_sub <-  pred_sub  %>%
        dplyr::filter(genotype %in% pull(dat_snp_all_unscaled, snp)) %>%
        dplyr::mutate(genotype = as.numeric(genotype)) %>%
        dplyr::mutate(!!snp := genotype)

   # Make PREDICTIONs
      
    # Full model  
      gam_preds <-  predict(gam_snp_int, 
                            newdata = pred_sub,
                            type = "response",
                            se = TRUE)
      pred_sub$pred <- gam_preds$fit
      pred_sub$se <- gam_preds$se.fit
      
   
   # Save which accession was used in prediction
       gam_snp_int["accession_pred"] <- as.character(pred_sub$accession[1])


    ## Write model summary to file
    # Save into dataframe 
    gam_mod_summary_df <- data.frame(
      
      # climate var
      climate_var = climate_var_dif,
      
      # Snp name
      snp = snp,
      
      # Accession used for prediction
      acc_pred = gam_snp_int$accession_pred,
      
      # Sample size
      n = gam_snp_int_summary$n,
      
      # Sample size of genotypes
      n_gen0 = gam_snp_int$n_gen0,
      n_gen1 = gam_snp_int$n_gen1,
      n_gen2 = gam_snp_int$n_gen2,

      # Deviance explained
      dev_explained = gam_snp_int_summary$dev.expl,
      
      # Difference in deviance
      dev_dif = gam_snp_int$dev_dif,
      
      # Rsquared adjusted
      rsq_adj = gam_snp_int_summary$r.sq,
      
      # Difference in r squared
      rsq_dif = gam_snp_int$rsq_dif,
      
      # # Use grep to look at the last number of the smooth summary table - 
      # # Should correspond to the genotype
       # # Returns NA if no match
     
     p_val_gen_int = ifelse(length(gam_snp_int_summary$s.table[grep(pattern = ":snp_",                                            rownames(gam_snp_int_summary$s.table)), "p-value"]) == 0,
                              yes = NA,
                              no = gam_snp_int_summary$s.table[grep(pattern = ":snp_",
                             rownames(gam_snp_int_summary$s.table)), "p-value"]),
     
     ## F value for interaction term
     f_val_gen_int = ifelse(length(gam_snp_int_summary$s.table[grep(pattern = ":snp_",                                            rownames(gam_snp_int_summary$s.table)), "F"]) == 0,
                            yes = NA,
                            no = gam_snp_int_summary$s.table[grep(pattern = ":snp_",
                                          rownames(gam_snp_int_summary$s.table)), "F"]),
     
     # P and f value for main term
     p_val_gen_main = ifelse(length(gam_snp_int_summary$p.table[grep(pattern = "snp_",                                            rownames(gam_snp_int_summary$p.table)), "Pr(>|t|)"]) == 0,
                            yes = NA,
                            no = gam_snp_int_summary$p.table[grep(pattern = "snp_",
                                      rownames(gam_snp_int_summary$p.table)), "Pr(>|t|)"]),
     
     ## T value for main term
     f_val_gen_main = ifelse(length(gam_snp_int_summary$p.table[grep(pattern = "snp_",                                            rownames(gam_snp_int_summary$p.table)), "t value"]) == 0,
                            yes = NA,
                            no = gam_snp_int_summary$p.table[grep(pattern = "snp_",
                 rownames(gam_snp_int_summary$p.table)), "t value"])
        )
    
   print(gam_mod_summary_df)
    
    
    # ## Write summary to file
    
    # If first time, write column names, if not, append to CSV
    if(x == 1){
      append = FALSE
    } else if(x > 1){
      append = TRUE
    }
    
      write_csv(gam_mod_summary_df, path = paste0("./model_summaries/gam_summaries_", task_id, ".csv"),
               append = append)
    

   ## Write predictions to file
      pred_sub_out <- pred_sub %>%
        dplyr::mutate(snp = snp) %>%
        dplyr::select(snp, tmax_dif_unscaled, genotype, pred, se)
      
      write_csv(pred_sub_out, path = paste0("./model_predictions/gam_predictions_", task_id, ".csv"),
               append = append)
    
    
    x = x + 1 # Increment loop
    
  # Timing of loop
    end_time <- Sys.time()
    cat("This loop took: ")
    print(end_time - start_time)
    cat("\n")   
    
  } # End SNP loop




