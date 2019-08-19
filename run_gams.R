# Read in command line arguments
	args<-commandArgs(trailingOnly=TRUE)

  print(args)

  task_id <- as.numeric(args[1])
  interval <- as.numeric(args[2])

# Load libraries
  library(mgcv, lib.loc = "/u/home/l/lukembro/R/x86_64-pc-linux-gnu-library/3.5") # Make sure to load latest version of mgcv
  library(tidyverse)
  # library(visreg) # Installed on home directory
  # library(patchwork) # Installed on home directory

# Load data
  load("../gam_cluster_2019-08-03.Rdata")
 
## Read in snps to skip
  skip_these_SNPs <- read_csv("../skip_these_SNPs_2019_08_05_acrossfams.csv")
  

# Initialize variables
  x = 1
 
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
     if(snp %in% skip_these_SNPs$snp){
       cat("Skipping SNP:", snp, " ... \n")
       next
     }
    
   # Make sure there's at least 2 levels
    if(length(table(dat_snp_all[, snp])) < 2){
      cat("Skipping because not at least 2 levels")
      next
    }

# Run gams
   	cat("Working on: ", snp, "... number: ", x, " ...\n" )

   ## Response variable as residual
   	fixed_effects_resids <- paste0(paste0("rgr_resids ~ ", snp, " + s(tmax_sum_dif, by = ", snp, ", bs=\"cr\")"))

   ## Impute missing data to mean allele frequency - assumes genotypes are scaled
    dat_snp_all[is.na(dat_snp_all[, snp]), snp] <- 0
     
  # Gam with interaction effect
    gam_snp_int = bam(
                      formula = formula(fixed_effects_resids),
                      data = dat_snp_all,
                      discrete = FALSE,
                      nthreads = 8,
                      method = "fREML")
                
    print(summary(gam_snp_int))
    
    gam_snp_int_summary <- summary(gam_snp_int)
    
    # Save sample size per genotype
    gen_table <-  table(dat_snp_all_unscaled[, snp])
    gam_snp_int$n_gen0 <- as.numeric(gen_table["0"])
    gam_snp_int$n_gen1 <- as.numeric(gen_table["1"])
    gam_snp_int$n_gen2 <- as.numeric(gen_table["2"])
    
    
  # Make predictions
    
    # Set up each climate variable individually
    pred <-  expand.grid(height_2014 = 0,
                         accession = "1",
                         section_block = "Block1_1",
                         locality = "FHL",
                         PC1_clim =0, PC2_clim = 0,
                         PC1_gen = 0, PC2_gen = 0,
                         tmax_sum_dif = seq(min(dat_snp_all$tmax_sum_dif), 
                                            max(dat_snp_all$tmax_sum_dif),
                                            length.out = 150))
    
    pred$tmax_sum_dif_unscaled <- back_transform(pred$tmax_sum_dif,
                                                    "tmax_sum_dif",
                                                    scaled_var_means_gbs_all,
                                                    scaled_var_sds_gbs_all)
    
    # When SNPs are continuous, make own genotypes
    gen_temp_all <- data.frame(accession = "1", 
                               genotype = c(0, 1, 2),             
                               genotype_scaled = c((0 - scaled_snps_means_all[snp]) /
                                                     scaled_snps_sds_all[snp],
                                                   (1 - scaled_snps_means_all[snp]) / scaled_snps_sds_all[snp],
                                                   (2 - scaled_snps_means_all[snp]) / scaled_snps_sds_all[snp]))
    
    # Join to prediction dataframe
    pred_sub <- left_join(pred, gen_temp_all)
    
    ## For SNPs as - make it so that column name is snp name
    pred_sub <-  pred_sub  %>%
      dplyr::filter(genotype %in% pull(dat_snp_all_unscaled, snp)) %>%
      dplyr::mutate(!!snp := genotype_scaled) %>% 
      dplyr::mutate(genotype = as.character(genotype))
    
    
    # Make prediction with full model
    gam_preds <-  predict(gam_snp_int, 
                          newdata = pred_sub,
                          type = "response",
                          se = TRUE)
    pred_sub$pred <- gam_preds$fit
    pred_sub$se   <- gam_preds$se.fit
    
    
    # Save which accession was used in prediction
    gam_snp_int["accession_pred"] <- as.character(pred_sub$accession[1])
    
    
    
## Format summary file
    # Save into dataframe 
    gam_mod_summary_df <- data.frame(
      
      # climate var
      climate_var = "tmax_sum_dif",
      
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
      
      # Rsquared adjusted
      rsq_adj = gam_snp_int_summary$r.sq,

      
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
                                                                   rownames(gam_snp_int_summary$s.table)), "F"])

    )
    
    gam_mod_summary_df
    
    
    


# Start Cross Validation loop ---------------------------------------------


for(fold in 1:length(training_folds)){
  cat(paste0("Working on fold: ", fold, "... \n"))
  
  ## Make training set
  ## Join residuals to training and testing data
  dat_snp_training <- dat_snp_all[training_folds[[fold]], ]
  dat_snp_training_unscaled <- dat_snp_all_unscaled[training_folds[[fold]], ]
  dim(dat_snp_training_unscaled)
  

  # Gam with interaction effect on training data
  gam_snp_int_training = bam(
                            formula = formula(fixed_effects_resids),
                            data = dat_snp_training,
                            discrete = FALSE,
                            nthreads = 8,
                            method = "fREML")
  
  print(summary(gam_snp_int_training)) 
  
  
  # Save summary
  gam_snp_int_training_summary <- summary(gam_snp_int_training)
  
  
  # Make prediction dataframe
  gen_temp_training <- data.frame(accession = "1", 
                                  genotype = c(0, 1, 2),             
                                  genotype_scaled = c((0 - scaled_snps_means_all[snp]) /
                                                        scaled_snps_sds_all[snp],
                                                      (1 - scaled_snps_means_all[snp]) / scaled_snps_sds_all[snp],
                                                      (2 - scaled_snps_means_all[snp]) / scaled_snps_sds_all[snp]))
  
  pred_sub_training <- left_join(pred, gen_temp_training)
  
  
  pred_sub_training <-  pred_sub_training  %>%
    dplyr::filter(genotype %in% pull(dat_snp_training_unscaled, snp)) %>%
    dplyr::mutate(!!snp := genotype_scaled) %>% 
    dplyr::mutate(genotype = as.character(genotype))
  
  
  
  # Make predictions

  # Training model  
  gam_preds_training <-  predict(gam_snp_int_training, 
                                 newdata = pred_sub_training,
                                 type = "response",
                                 se = TRUE)
  
  pred_sub_training[, paste0("pred_fold", fold)] <- gam_preds_training$fit

  # Join training data to full model predictions - should result in adding predictions for this fold as a new column
  pred_sub <-  dplyr::left_join(pred_sub,
                                dplyr::select(pred_sub_training, genotype,
                                              tmax_sum_dif, 
                                              paste0("pred_fold", fold))) 
  
  
  
  ## Extract test statistics and p value for training model
  p_val_gen_int_training = ifelse(length(gam_snp_int_training_summary$s.table[grep(pattern = ":snp_",
                                 rownames(gam_snp_int_training_summary$s.table)), "p-value"]) == 0,
                                  yes = NA,
                                  no = gam_snp_int_training_summary$s.table[grep(pattern = ":snp_",
                                         rownames(gam_snp_int_training_summary$s.table)), "p-value"])
  
  ## F value for interaction term
  f_val_gen_int_training = ifelse(length(gam_snp_int_training_summary$s.table[grep(pattern = ":snp_", 
                                           rownames(gam_snp_int_training_summary$s.table)), "F"]) == 0,
                                  yes = NA,
                                  no = gam_snp_int_training_summary$s.table[grep(pattern = ":snp_",
                                                     rownames(gam_snp_int_training_summary$s.table)), "F"])
  
  ## Add to summary df
  gam_mod_summary_df[, paste0("p_val_gen_int_fold", fold)] <- p_val_gen_int_training
  gam_mod_summary_df[, paste0("f_val_gen_int_fold", fold)] <- f_val_gen_int_training
  
} ## End fold loop   
    
    
  

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
        # Get rid of extra columns
        dplyr::select(-height_2014, -accession, -genotype_scaled, -section_block,
                      -locality, -PC1_clim, -PC2_clim, -PC1_gen, -PC2_gen, -!!snp) %>%
        # Put columns in right order
        dplyr::select(snp, genotype, tmax_sum_dif, everything())
      
      write_csv(pred_sub_out, path = paste0("./model_predictions/gam_predictions_", task_id, ".csv"),
               append = append)
    
    
    x = x + 1
    
    # Timing of loop
    end_time <- Sys.time()
    cat("This loop took: ")
    print(end_time - start_time)
    cat("\n")   
    
  } # End SNP loop




