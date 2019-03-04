# Read in command line arguments
	args<-commandArgs(trailingOnly=TRUE)

  print(args)

	task_id <- as.numeric(args[1])
  interval <- as.numeric(args[2])
  climate_var_dif <- as.character(args[3])

# Load libraries
  library(mgcv, lib.loc = "/u/home/l/lukembro/R/x86_64-pc-linux-gnu-library/3.5") # Make sure to load latest version of mgcv
  library(tidyverse)
  # library(visreg) # Installed on home directory
  # library(patchwork) # Installed on home directory

# Load data
  load("../gam_cluster_2019-02-26.Rdata")
 
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
    random_effects <-  '+ s(accession, bs = "re") + s(locality, bs = "re")'
  

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
    if(length(table(dat_snp_training_unscaled[, snp])) < 2){
      cat("Skipping because not at least 2 levels")
      next
    }

# Run gams
   	cat("Working on: ", snp, "... number: ", x, " ...\n" )

  # For SNPs as continuous
   #   fixed_effects_int <- paste0("rgr ~ section_block + ", snp, " + s(height_2014, bs =\"cr\") + s(", climate_var_dif, " , bs = \"cr\") + s(", climate_var_dif,", by = ", snp, ", bs=\"cr\") + s(PC1_clim, bs =\"cr\") + s(PC2_clim, bs =\"cr\") + s(tmax_sum_dif, by = PC1_clim, bs =\"cr\") + s(tmax_sum_dif, by = PC2_clim, bs =\"cr\")")
   # 
   #  ## Add gen pcs and interactions
   # fixed_effects_int <- paste0(fixed_effects_int, paste0(paste0("+ s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC1_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC2_gen, bs=\"cr\")")))
   	
   	
   ## Testing out residual gams
   	fixed_effects_resids <- paste0(paste0("rgr_resids ~ ", snp, " + s(tmax_sum_dif, by = ", snp, ", bs=\"cr\")"))
   	

   ## Randomize RGR data
    # dat_snp_training_unscaled$rgr <- sample(dat_snp_training_unscaled$rgr)

   ## Randomize Genotype
   #  dat_snp_training[, snp] <- sample(pull(dat_snp_training, snp))
    #  dat_snp_training_unscaled[, snp] <- sample(x = c(0, 1 , 2), 
    #                                    size = nrow(dat_snp_training_unscaled), 
    #                                    replace = TRUE)
   
   ## Impute missing data to mean allele frequency - assumes genotypes are scaled
    dat_snp_training[is.na(dat_snp_training[, snp]), snp] <- 0
     
  # Gam with interaction effect
    gam_snp_int = bam(
                 # formula = make_formula(fixed_effects_int, random_effects),
                   formula = formula(fixed_effects_resids),
                  data = dat_snp_training,
                  discrete = TRUE,
                  nthreads = 8,
                  method = "fREML")
                #  family = "tw")

    print(summary(gam_snp_int))


 ## Calculate decrease in deviance with SNP data  
    
    # # Fit model with no SNP data but with PCs
    #   form_no_snp <- formula(paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"cr\") + s(", climate_var_dif,", bs=\"cr\") + s(locality, bs = \"re\") + s(accession, bs = \"re\") + s(tmax_sum_dif, by = PC1_clim, bs = \"cr\") + s(tmax_sum_dif, by = PC2_clim, bs = \"cr\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC1_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC2_gen, bs=\"cr\")")))
    # 
    #   gam_no_snp = bam(formula = form_no_snp,
    #                     data = dat_snp_training,
    #                     discrete = TRUE,
    #                     nthreads = 8,
    #                     method = "fREML",
    #                     family = "tw",
    #                     control = list(trace = FALSE))
    #   
    # #  summary(gam_no_snp)
    #   
   
  # Save deviance differences and Rsquared into model
    gam_snp_int_summary <- summary(gam_snp_int)
  #  gam_no_snp_summary  <- summary(gam_no_snp)
    
 #   gam_snp_int$dev_dif <- gam_snp_int_summary$dev - gam_no_snp_summary$dev
  #  gam_snp_int$rsq_dif <- gam_snp_int_summary$r.sq - gam_no_snp_summary$r.sq
    
    gam_snp_int$dev_dif <- NA
    gam_snp_int$rsq_dif <- NA
    
   # Save sample size per genotype
     gen_table <-  table(dat_snp_training_unscaled[, snp])
     gam_snp_int$n_gen0 <- as.numeric(gen_table["0"])
     gam_snp_int$n_gen1 <- as.numeric(gen_table["1"])
     gam_snp_int$n_gen2 <- as.numeric(gen_table["2"])
    

  # Make predictions

   # When SNPs are continuous, make own genotypes
      gen_temp <- data.frame(accession = "1", 
                 genotype = c(0, 1, 2),             
               	 genotype_scaled = c((0 - scaled_snps_means_training[snp]) /
               	                       scaled_snps_sds_training[snp],
                              (1 - scaled_snps_means_training[snp]) / scaled_snps_sds_training[snp],
                              (2 - scaled_snps_means_training[snp]) / scaled_snps_sds_training[snp]))
                 
      
      pred_sub <- left_join(pred, gen_temp)
        
        
    # Make sure genotypes are represented and rename column      
    # pred_sub <-  pred_sub  %>%
    #   dplyr::filter(genotype_scaled %in% pull(dat_snp_training, snp)) %>% 
    #   dplyr::mutate(!!snp := genotype_scaled)

  ## For SNPs as factors
      # pred_sub <-  pred_sub  %>%
      #   dplyr::filter(genotype %in% pull(dat_snp_training_unscaled, snp)) %>%
      #   dplyr::mutate(!!snp := as.character(genotype)) %>% # Change to character for factoring
      #   dplyr::mutate(genotype = as.character(genotype))
      
  ## For SNPs as continuous
      pred_sub <-  pred_sub  %>%
        dplyr::filter(genotype %in% pull(dat_snp_training_unscaled, snp)) %>%
        dplyr::mutate(!!snp := genotype_scaled) %>% 
        dplyr::mutate(genotype = as.character(genotype))

	# Change accession number if genetic data is missing so we don't get errors in prediction
        # while(!(pred_sub$accession[1] %in%  gam_snp_int$model$accession)){
        #   pred_sub$accession <- as.character(as.numeric(pred_sub$accession) + 1)
        # }

  ## Visreg plots
  #  pdf(paste0("./model_plots/", snp,".pdf"), width = 8, height = 5)
#
  # # With partial residuals
  #  visreg(gam_snp_int, xvar = climate_var_dif, by = snp,
  #         overlay = TRUE, partial = TRUE, 
  #         breaks = c(unique(pred_sub$genotype_scaled)),
  #         xtrans = function(x) {(x *  scaled_var_sds_gbs_only[climate_var_dif]) + scaled_var_means_gbs_only[climate_var_dif]},
  #         ylab = "5 year height (cm)")
  #  
  #  # Without partial residuals
   # visreg(gam_snp_int, xvar = climate_var_dif, by = snp,
   #        overlay = TRUE, partial = FALSE, rug = FALSE,
   #        breaks = c(unique(pred_sub$genotype_scaled)),
   #        xtrans = function(x) {(x *  scaled_var_sds_gbs_only[climate_var_dif]) + scaled_var_means_gbs_only[climate_var_dif]},
   #        ylab = "5 year height (cm)")
  #  
  #  dev.off()      


   # Make PREDICTIONs
      
    # Full model  
      gam_preds <-  predict(gam_snp_int, 
                            newdata = pred_sub,
                            type = "response",
                            se = TRUE)
      pred_sub$pred <- gam_preds$fit
      pred_sub$se <- gam_preds$se.fit

     
### Training and testing validation
      
      
    # # Select just columns that we will use to speed it up
    # dat_snp_training_reduced <- dplyr::select(dat_snp_training, rgr, section_block, snp,
    #                               height_2014, tmax_sum_dif,
    #                               PC1_clim, PC2_clim, PC1_gen, PC2_gen, locality, accession)
    # 
    # for(fold in 1:length(folds)){
    # 
    #     dat_snp_training_train <- dat_snp_training_reduced[-folds[[fold]], ]
    #     dat_snp_training_train <-  dat_snp_training_train[!is.na(pull(dat_snp_training_train, snp)), ]
    # 
    #    # dim(dat_snp_training_train)
    # 
    #     dat_snp_training_test <- dat_snp_training_reduced[folds[[fold]], ]
    #     dat_snp_training_test <-  dat_snp_training_test[!is.na(pull(dat_snp_training_test, snp)), ]
    # 
    #     # Make sure all accessions are represented to avoid prediction errors
    #     dat_snp_training_test <- dat_snp_training_test[dat_snp_training_test$accession %in% dat_snp_training_train$accession, ]
    # 
    #     # dim(dat_snp_training_test)
    # 
    #   # Training model
    #    gam_snp_int_train = bam(formula = make_formula(fixed_effects_int, random_effects),
    #                             data = dat_snp_training_train,
    #                             discrete = TRUE,
    #                             nthreads = 8,
    #                             method = "fREML",
    #                             family = "tw")
    # 
    #     # print(summary(gam_snp_int_train))
    # 
    #   # Calculate predicted growth rates of test set based on training set
    #     test_rgr = predict(gam_snp_int_train,
    #                    newdata = dat_snp_training_test,
    #                    type = "response")
    # 
    #     # Calculate correlation separately for each genotype
    #     dat_snp_training_test$rgr_predicted <- test_rgr
    # 
    #    cv_cors_temp <-  dat_snp_training_test %>%
    #       dplyr::select(snp, rgr, rgr_predicted) %>%
    #       group_by(.dots = snp) %>% # Group by genotype
    #      summarise(cv_cor = cor(rgr, rgr_predicted))
    # 
    #    if(fold == 1){
    #      cv_cors = cv_cors_temp
    #    } else {
    #      cv_cors <- rbind(cv_cors, cv_cors_temp)
    #    }
    # 
    #   } # End cross validation loop
    # 
    # ## Average across cross validation folds
    # cv_cors_avg <- cv_cors %>%
    #           group_by(.dots = snp) %>%
    #           summarise(cv_cor = mean(cv_cor))
    # 
    # cv_cors_avg[, 1] <- back_transform(as.data.frame(cv_cors_avg[, 1]), var = snp,
    #                                    means = scaled_snps_means, sds = scaled_snps_sds)
    # 
    # cv_cor_gen_0 <- as.numeric(cv_cors_avg[pull(cv_cors_avg, snp) == 0, 2])
    # cv_cor_gen_1 <- as.numeric(cv_cors_avg[pull(cv_cors_avg, snp) == 1, 2])
    # cv_cor_gen_2 <- as.numeric(cv_cors_avg[pull(cv_cors_avg, snp) == 2, 2])
    #  
      
    # If not running CV
       cv_cor_gen_0 = NA
       cv_cor_gen_1 = NA
       cv_cor_gen_2 = NA
   
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
      
      # Converged?
      converged = gam_snp_int$mgcv.conv,
      
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
                 rownames(gam_snp_int_summary$p.table)), "t value"]),
     
     
     
     # Cross validation correlation for each genotype
     cv_cor_gen_0 = cv_cor_gen_0,
     cv_cor_gen_1 = cv_cor_gen_1,
     cv_cor_gen_2 = cv_cor_gen_2
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
        dplyr::select(snp, tmax_sum_dif_unscaled, genotype, pred, se)
      
      write_csv(pred_sub_out, path = paste0("./model_predictions/gam_predictions_", task_id, ".csv"),
               append = append)
    
    
    x = x + 1
    
    
  # Save individual .Rdata file for each model
   # cat("Saving gam model to file... \n")

   # gam_name <- paste0("gam_", climate_var_dif, "_", snp)
   # assign(gam_name, gam_snp_int) # Assign model the name

    # save(list = gam_name, file = paste0("./gam_mods_data/", gam_name, ".Rdata"))
    
    # Timing of loop
    end_time <- Sys.time()
    cat("This loop took: ")
    print(end_time - start_time)
    cat("\n")   
    
  } # End SNP loop




