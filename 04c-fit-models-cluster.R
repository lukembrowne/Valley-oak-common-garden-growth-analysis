

# Run model setup ---------------------------------------------------------

  source("./04b-setup-models-cluster.R")


# Gam code from cluster for testing -------------------------------------------

# Initialize variables
  gam_list <- list()
  x = 1
  
  climate_var_dif = "tmax_sum_dif"
  
  ## Convert fixed and random parts of formula to actual formula object
  # No idea why we need to call formula twice, but it works
  make_formula <- function(fixed, random){
    return(formula(formula(enquote(paste0(fixed, random)))))
  }
  
   # Formula for random effects
  random_effects <-  '+ s(accession, bs = "re")'
  
  
  task_id = 5; interval = 5
  save_outlier_preds <- FALSE
  outlier_preds <- NULL

# # # Loop through snps by index based on task id
  for(snp_index in task_id:(task_id + interval - 1)){

    # To avoid going past number of snps
    if(snp_index > length(snp_col_names)){
      next
    }

    # Choose snp
    snp <- snp_col_names[snp_index]

 # # LOOPING THROUGH TOP SNPS
 # for(snp in c(top_snps_long$snp, bottom_snps_long$snp, mid_snps_long$snp)){
 #   save_outlier_preds <- TRUE # Save output of the top SNPs for plotting together
 # 
    
    
      # For timing loops
      start_time <- Sys.time()
    

    # Make sure there's at least 2 levels
    if(length(table(dat_snp_unscaled[, snp])) < 2){
      cat("Skipping because not at least 2 levels")
      next
    }
    
    # Convert to factor
   # dat_snp[, snp] <- as.factor(pull(dat_snp, snp))
    dat_snp_unscaled[, snp] <- as.factor(pull(dat_snp_unscaled, snp))
    
     # Formula for fixed effects
    # Stack overflow on adding m = 1 to by= smooths - https://stats.stackexchange.com/questions/32730/how-to-include-an-interaction-term-in-gam

    
# Run gams  
    cat("Working on: ", snp, "... number: ", x, " ...\n" )    
    
  # For SNPS as continuous  
   # fixed_effects_int <- paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"cr\") + te(", snp, ", bs=\"cr\", k = ", k, ") +  te(", climate_var_dif,", bs=\"cr\") + ti(", climate_var_dif,", ",snp,", bs=\"cr\", k =", k," )"))
    
    fixed_effects_int <- paste0(paste0("rgr ~ section_block + ", snp, "+ s(height_2014, bs=\"cr\") +  s(", climate_var_dif,", bs=\"cr\") + s(", climate_var_dif,", by = ", snp, ", bs=\"cr\")"))
    
    ## Add gen pcs and interactions
    fixed_effects_int <- paste0(fixed_effects_int, paste0(paste0("+ s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC1_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC2_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC3_gen, bs=\"cr\")")))
    
    ## Add MEMs
  #  fixed_effects_int <- paste0(fixed_effects_int,  "+ s(mem1, bs =\"cr\") + s(mem2, bs =\"cr\") + s(mem3, bs =\"cr\") + s(mem4, bs =\"cr\")")
    
  ## Randomize RGR data
   # dat_snp_unscaled$rgr <- sample(dat_snp_unscaled$rgr)

   # Gam with interaction effect
    gam_snp_int = bam(formula = make_formula(fixed_effects_int, random_effects),
                  data = dat_snp_unscaled[!is.na(pull(dat_snp_unscaled, snp)), ],
                 discrete = TRUE, # Discrete makes the right curve up and leads to low p values!!
                  nthreads = 8,
                method = "fREML",
             #   method = "ML",
                  family = "tw",
                  control = list(trace = FALSE))

    print(summary(gam_snp_int))
    
  
  ## Calculate decrease in deviance with SNP data  
    
    # Fit model with no SNP data but with PCs
      form_no_snp <- formula(paste0("rgr ~ section_block + s(height_2014, bs =\"cr\") + s(", climate_var_dif, " , bs = \"cr\") + s(accession, bs = \"re\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC1_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC2_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC3_gen, bs=\"cr\")"))

      gam_no_snp = bam(formula = form_no_snp,
                        data = dat_snp_unscaled[!is.na(pull(dat_snp_unscaled, snp)), ],
                        discrete = TRUE,
                        nthreads = 8,
                        method = "fREML",
                        #   method = "ML",
                        family = "tw",
                        control = list(trace = FALSE))

    #  summary(gam_no_snp)

    # Adding random variable
    form_random <-formula(paste0(paste0("rgr ~ section_block + random + s(height_2014, bs=\"cr\") +  s(", climate_var_dif,", bs=\"cr\") + + s(accession, bs = \"re\") + s(", climate_var_dif,", by = random , bs=\"cr\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC1_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC2_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC3_gen, bs=\"cr\")")))

    gam_random = bam(formula = form_random,
                     data = dat_snp_unscaled[!is.na(pull(dat_snp_unscaled, snp)), ],
                     discrete = TRUE,
                     nthreads = 8,
                     method = "fREML",
                     #   method = "ML",
                     family = "tw",
                     control = list(trace = FALSE))

   # summary(gam_random)

  # Save deviance differences and Rsquared into model
    gam_snp_int_summary <- summary(gam_snp_int)
    gam_no_snp_summary  <- summary(gam_no_snp)
    gam_random_summary  <- summary(gam_random)

    gam_snp_int$dev_dif <- gam_snp_int_summary$dev - gam_no_snp_summary$dev
    gam_snp_int$dev_dif_random <-  gam_random_summary$dev - gam_no_snp_summary$dev

    gam_snp_int$rsq_dif <- gam_snp_int_summary$r.sq - gam_no_snp_summary$r.sq
    gam_snp_int$rsq_dif_random <- gam_random_summary$r.sq - gam_no_snp_summary$r.sq
    
    
   ## Save sample size per genotype
      gen_table <-  table(dat_snp_unscaled[, snp])
      
      gam_snp_int$n_gen0 <- as.numeric(gen_table["0"])
      gam_snp_int$n_gen1 <- as.numeric(gen_table["1"])
      gam_snp_int$n_gen2 <- as.numeric(gen_table["2"])
      
      
  # Visualization plots
    
   # For snps as factors
    # pdf(paste0("./output/model_visualizations/", snp,".pdf"), width = 8, height = 5)
    # visreg(gam_snp_int, xvar = climate_var_dif, overlay = TRUE, main = snp,
    #        by = snp, scale = "response")
    # dev.off()
    
    # For SNPS as continuous
        # pdf(paste0("./output/model_visualizations/", snp,".pdf"), width = 8, height = 5)
        # visreg2d(gam_snp_int, 
        #          xvar = climate_var_dif, yvar = snp, 
        #          scale = "response", 
        #          plot.type = "persp", theta = 35)
        # 
        # visreg(gam_snp_int, xvar = climate_var_dif, by = snp, breaks = gen_temp$genotype_scaled, scale = "response")
        # 
        # 
        # dev.off()
    
    # Make predictions
        
    # When SNPs are continuous, make own genotypes
      gen_temp <- data.frame(accession = "1", 
                  genotype = c(0, 1, 2),             
                  genotype_scaled = c((0 - scaled_snps_means[snp]) / scaled_snps_sds[snp],
                              (1 - scaled_snps_means[snp]) / scaled_snps_sds[snp],
                              (2 - scaled_snps_means[snp]) / scaled_snps_sds[snp]))
                 
      
      pred_sub <- left_join(pred, gen_temp)
        
        
    # Make sure genotypes are represented and rename column      
      # pred_sub <-  pred_sub  %>%
      #   dplyr::filter(genotype_scaled %in% pull(dat_snp, snp)) %>%
      #   dplyr::mutate(!!snp := genotype_scaled)
    
      ## For SNPs as factors
      pred_sub <-  pred_sub  %>%
        dplyr::filter(genotype %in% pull(dat_snp_unscaled, snp)) %>%
        dplyr::mutate(!!snp := as.character(genotype)) %>% # Change to character for factoring
        dplyr::mutate(genotype = as.character(genotype))
      
      # Change accession number if genetic data is missing so we don't get errors in prediction
        while(!(pred_sub$accession[1] %in%  gam_snp_int$model$accession)){
          pred_sub$accession <- as.character(as.numeric(pred_sub$accession) + 1)
        }
      
    ## Remove extra columns
      pred_sub <- pred_sub %>%
        dplyr::select(-tmin_dif_unscaled, -tmin_dif,
                      -tave_dif_unscaled, -tave_dif,
                      -tave_lgm_dif_unscaled, -tave_lgm_dif,
                      -tmax_sum_lgm_dif_unscaled, -tmax_sum_lgm_dif,
                      -tmin_winter_dif_unscaled, -tmin_winter_dif,          
                     -tmin_winter_lgm_dif_unscaled, -tmin_winter_lgm_dif,         
                     -DD5_dif_unscaled, -DD5_dif,                    
                     -DD5_lgm_dif_unscaled, -DD5_lgm_dif,             
                     -random_dif_unscaled, -random_dif)
      
      

    # Make PREDICTIONs
      
    # Full model  
      gam_preds <-  predict(gam_snp_int, 
                            newdata = pred_sub,
                            type = "response",
                            se = TRUE)
      pred_sub$pred <- gam_preds$fit
      pred_sub$se <- gam_preds$se.fit
      
    
      
      
### Training and testing validation
      
    for(fold in 1:length(folds)){

        dat_snp_train <- dat_snp_unscaled[-folds[[fold]], ]      
        dat_snp_train <-  dat_snp_train[!is.na(pull(dat_snp_train, snp)), ]
        
        dim(dat_snp_train)
        
        dat_snp_test <- dat_snp_unscaled[folds[[fold]], ]
        dat_snp_test <-  dat_snp_test[!is.na(pull(dat_snp_test, snp)), ]
        dim(dat_snp_test)
        
      # Training model  
       gam_snp_int_train = bam(formula = make_formula(fixed_effects_int, random_effects),
                                data =dat_snp_train,
                                discrete = TRUE, 
                                nthreads = 8,
                                method = "fREML",
                                family = "tw")
        
        # print(summary(gam_snp_int_train))
      # Testing model  
        gam_snp_int_test = bam(formula = make_formula(fixed_effects_int, random_effects),
                               data = dat_snp_test,
                               discrete = TRUE, 
                               nthreads = 8,
                               method = "fREML",
                               family = "tw")
        
        # print(summary(gam_snp_int_test))
        
        
        # Save summaries
        # gam_snp_int_train_summary <- summary(gam_snp_int_train)
        # gam_snp_int_test_summary <- summary(gam_snp_int_test)
        
      
      
      ## Making predictions
 
        # Change accession number if genetic data is missing so we don't get errors in prediction
        while(!(pred_sub$accession[1] %in%  gam_snp_int_train$model$accession) | 
              !(pred_sub$accession[1] %in%  gam_snp_int_test$model$accession)){
          pred_sub$accession <- as.character(as.numeric(pred_sub$accession) + 1)
        }
        
        
        # Filter down to exclude gneotypes that are not in model to avoid error in prediction
        pred_sub_train <- pred_sub %>%
                          dplyr::filter(genotype %in% pull(dat_snp_train, snp))
        table(pred_sub_train$genotype)
        
        gam_preds_train <-  predict(gam_snp_int_train, 
                                    newdata = pred_sub_train,
                                    type = "response",
                                    se = TRUE)
        
        pred_sub_train[paste0("pred_train_", fold)] <- gam_preds_train$fit
        pred_sub_train[paste0("se_train_", fold)] <- gam_preds_train$se.fit
        
        # Join back to main prediction df
        pred_sub <- left_join(pred_sub, 
                             dplyr::select(pred_sub_train, 
                                           genotype, tmax_sum_dif, 
                                           paste0("pred_train_", fold),
                                           paste0("se_train_", fold)))
        
        
      ## Testing model
       
       # Filter down to exclude gneotypes that are not in model to avoid error in prediction
       pred_sub_test <- pred_sub %>%
         dplyr::filter(genotype %in% pull(dat_snp_test, snp))
       table(pred_sub_test$genotype)
       
       
        gam_preds_test <-  predict(gam_snp_int_test, 
                                   newdata = pred_sub_test,
                                   type = "response",
                                   se = TRUE)
        
        pred_sub_test[paste0("pred_test_", fold)] <- gam_preds_test$fit
        pred_sub_test[paste0("se_test_", fold)] <- gam_preds_test$se.fit
        
        # Join back to main prediction df
        pred_sub <- left_join(pred_sub, 
                              dplyr::select(pred_sub_test, 
                                            genotype, tmax_sum_dif, 
                                            paste0("pred_test_", fold),
                                            paste0("se_test_", fold)))
      }
    
      # Factor genotypes
      pred_sub$genotype <- factor(pred_sub$genotype)
    
    # Backtransform climate variable
    # dat_snp_trans <- dat_snp_unscaled %>%
    #   dplyr::select(climate_var_dif, snp, rgr) %>%
    #   # Back transform X axis
    #   mutate(!!climate_var_dif := back_transform(x = get(climate_var_dif),
    #                                              var = climate_var_dif,
    #                                              means = scaled_var_means_gbs_only,
    #                                              sds = scaled_var_sds_gbs_only)) %>%
    #   # Filter down to just 0,1,2 genotypes
    #   dplyr::filter(get(snp) %in% c(0, 1, 2))
    
    
    # # Climate transfer functions by genotype
    #   p1 <-  
    #     pred_sub %>%
    #   #  dplyr::filter(genotype !=  "1") %>%
    #     ggplot(.) +
    #     geom_ribbon(aes(tmax_sum_dif_unscaled, ymin = pred - 1.96 * se, 
    #                     ymax = pred + 1.96 * se, fill = factor(genotype)),
    #                 alpha = 0.25, show.legend = FALSE)  +
    #     geom_vline(aes(xintercept = 0), lty = 2) +
    #     geom_vline(xintercept = 1.1) + 
    #     geom_vline(xintercept = 4.8, lwd = 1.5) +
    #     geom_line(aes(x = tmax_sum_dif_unscaled, y = pred, col = factor(genotype)), lwd = 1.5) +
    #     
    #     labs(col = "Genotype") +
    #     xlim(-2.5, 7.5) +
    #    scale_color_discrete(labels = c("Homozygous \n reference",
    #                                    "Heterozygous", "Homozygous\n non-reference")) +
    #     ylab("Relative growth rate") +
    #     xlab("Tmax transfer distance") +
    #     ggtitle(snp) +
    #     theme_bw(15) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    #                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    # 
    # # Density plots of genotypes
    #   p2 <- ggplot(dat_snp_trans,
    #                aes(x = get(climate_var_dif), fill = factor(get(snp)))) +
    #     geom_histogram(binwidth = .5, col = "grey10") +
    #     geom_vline(aes(xintercept = 0), lty = 2) +
    #     labs(fill = "Genotype") +
    #     xlab(climate_var_dif) +
    #     ggtitle(snp) +
    #     theme_bw(15)
    # 
    # # Plot to PDF
    #   pdf(paste0("./output/model_visualizations/", snp,"_gg.pdf"), width = 15, height = 5)
    #    print(p1 + p2)
    #   dev.off()
    #   
    
  # ## Save top and bottom predictions
  #     if(save_outlier_preds){
  #       
  #       type <- NA
  #       gen <- NA
  #       
  #       ## Figure out if top or bottom snp
  #       if(snp %in% top_snps_long$snp){
  #         type <- "top"
  #         gen <- as.numeric(top_snps_long$genotype[top_snps_long$snp == snp])
  #       } else if(snp %in% bottom_snps_long$snp){
  #         type <- "bottom"
  #         gen <- as.numeric(bottom_snps_long$genotype[bottom_snps_long$snp == snp])
  #       } else if(snp %in% mid_snps_long$snp){
  #         type <- "mid"
  #         gen <- as.numeric(mid_snps_long$genotype[mid_snps_long$snp == snp])
  #       }
  #       
  #      #  error check
  #       if(is.na(type) | is.na(gen)){
  #         stop("type or genotype is NA and it shouldn't be")
  #       }
  # 
  #      outlier_preds_temp <-  pred_sub %>% 
  #         dplyr::filter(get(snp) == gen) %>%
  #         mutate_if(is.factor, as.character) %>%
  #         mutate(snp = snp) %>%
  #         mutate(type = type) %>%
  #         dplyr::select(snp, type, tmax_sum_dif_unscaled, pred, se)
  #       
  #       if(!exists("outlier_preds")){
  #         outlier_preds <- outlier_preds_temp
  #       } else {
  #         outlier_preds <- rbind(outlier_preds, outlier_preds_temp)
  #       }
  # 
  #     } # End save outlier_preds
      
      

    
  # Make prediction for XXXX degree increase in temperature
      for(genotype in c(0,1,2)){
        
        if(!genotype %in% as.numeric(as.character(pred_sub$genotype))){
          gam_snp_int[paste0("height_4.8_gen_", genotype)] <- NA
          gam_snp_int[paste0("height_0_gen_", genotype)] <- NA
          gam_snp_int[paste0("height_change_gen_", genotype)] <- NA
          next
        }
        
        pred_sub_temp <- pred_sub[pred_sub$genotype == genotype, ]
        
        # Calculate predicted height
        height_4.8 <- pred_sub_temp$pred[pred_sub_temp$tmax_sum_dif_unscaled == 4.8][1]
        height_0 <- pred_sub_temp$pred[pred_sub_temp$tmax_sum_dif_unscaled == 0][1]
        
        # % Height change at 4.8
        height_change <- (height_4.8 - height_0) / height_0 * 100
           
        # Save values into model list
        gam_snp_int[paste0("height_4.8_gen_", genotype)] <- height_4.8
        gam_snp_int[paste0("height_0_gen_", genotype)] <- height_0
        gam_snp_int[paste0("height_change_gen_", genotype)] <- height_change
        
        # Save which accession was used in prediction
        gam_snp_int["accession_pred"] <- as.character(pred_sub$accession[1])
        
      }

    
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
      dev_dif_random = gam_snp_int$dev_dif_random,
      
      # Rsquared adjusted
      rsq_adj = gam_snp_int_summary$r.sq,
      
      # Difference in r squared
      rsq_dif = gam_snp_int$rsq_dif,
      rsq_dif_random = gam_snp_int$rsq_dif_random,
      
      # Predictions of change in height
      height_change_gen_0 = gam_snp_int$height_change_gen_0,
      height_change_gen_1 = gam_snp_int$height_change_gen_1,
      height_change_gen_2 = gam_snp_int$height_change_gen_2,
      
      # Predictions of absolute height
      height_4.8_gen_0 = gam_snp_int$height_4.8_gen_0,
      height_4.8_gen_1 = gam_snp_int$height_4.8_gen_1,
      height_4.8_gen_2 = gam_snp_int$height_4.8_gen_2,
      
      # Predictions of absolute height
      height_0_gen_0 = gam_snp_int$height_0_gen_0,
      height_0_gen_1 = gam_snp_int$height_0_gen_1,
      height_0_gen_2 = gam_snp_int$height_0_gen_2,
      
      
      # # Use grep to look at the last number of the smooth summary table - 
      # # Should correspond to the genotype
      # # Returns NA if no match
      p_val_gen_0 = ifelse(length(gam_snp_int_summary$s.table[grep(pattern = "0$",
                               rownames(gam_snp_int_summary$s.table)), "p-value"]) == 0,
                              yes = NA,
                              no = gam_snp_int_summary$s.table[grep(pattern = "0$",
                               rownames(gam_snp_int_summary$s.table)), "p-value"]),
      
      p_val_gen_1 = ifelse(length(gam_snp_int_summary$s.table[grep(pattern = "1$",
                          rownames(gam_snp_int_summary$s.table)), "p-value"]) == 0,
                           yes = NA,
                           no = gam_snp_int_summary$s.table[grep(pattern = "1$",
                           rownames(gam_snp_int_summary$s.table)), "p-value"]),
      
      p_val_gen_2 = ifelse(length(gam_snp_int_summary$s.table[grep(pattern = "2$",
                           rownames(gam_snp_int_summary$s.table)), "p-value"]) == 0,
                           yes = NA,
                           no = gam_snp_int_summary$s.table[grep(pattern = "2$",
                            rownames(gam_snp_int_summary$s.table)), "p-value"])
      
      
      # ## Training model
      # p_val_gen_0_train = ifelse(length(gam_snp_int_train_summary$s.table[grep(pattern = "0$",
      #                        rownames(gam_snp_int_train_summary$s.table)), "p-value"]) == 0,
      #                      yes = NA,
      #                      no = gam_snp_int_train_summary$s.table[grep(pattern = "0$",
      #                         rownames(gam_snp_int_train_summary$s.table)), "p-value"]),
      # 
      # p_val_gen_1_train = ifelse(length(gam_snp_int_train_summary$s.table[grep(pattern = "1$",
      #                       rownames(gam_snp_int_train_summary$s.table)), "p-value"]) == 0,
      #                      yes = NA,
      #                      no = gam_snp_int_train_summary$s.table[grep(pattern = "1$",
      #                        rownames(gam_snp_int_train_summary$s.table)), "p-value"]),
      # 
      # p_val_gen_2_train = ifelse(length(gam_snp_int_train_summary$s.table[grep(pattern = "2$",
      #                         rownames(gam_snp_int_train_summary$s.table)), "p-value"]) == 0,
      #                      yes = NA,
      #                      no = gam_snp_int_train_summary$s.table[grep(pattern = "2$",
      #                       rownames(gam_snp_int_train_summary$s.table)), "p-value"]),
      # 
      # ## Testing model
      # p_val_gen_0_test = ifelse(length(gam_snp_int_test_summary$s.table[grep(pattern = "0$",
      #                       rownames(gam_snp_int_test_summary$s.table)), "p-value"]) == 0,
      #                      yes = NA,
      #                      no = gam_snp_int_test_summary$s.table[grep(pattern = "0$",
      #                         rownames(gam_snp_int_test_summary$s.table)), "p-value"]),
      # 
      # p_val_gen_1_test = ifelse(length(gam_snp_int_test_summary$s.table[grep(pattern = "1$",
      #                        rownames(gam_snp_int_test_summary$s.table)), "p-value"]) == 0,
      #                      yes = NA,
      #                      no = gam_snp_int_test_summary$s.table[grep(pattern = "1$",
      #                         rownames(gam_snp_int_test_summary$s.table)), "p-value"]),
      # 
      # p_val_gen_2_test = ifelse(length(gam_snp_int_test_summary$s.table[grep(pattern = "2$",
      #                         rownames(gam_snp_int_test_summary$s.table)), "p-value"]) == 0,
      #                      yes = NA,
      #                      no = gam_snp_int_test_summary$s.table[grep(pattern = "2$",
      #                        rownames(gam_snp_int_test_summary$s.table)), "p-value"])
      # 
      
      
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
    

    # ## Write predictions to file
      pred_sub_out <- pred_sub %>%
        dplyr::mutate(snp = snp) %>%
        dplyr::select(snp, tmax_sum_dif_unscaled, genotype, pred, se,
                      contains("train"),
                      contains("test"))

      write_csv(pred_sub_out, path = paste0("./model_predictions/gam_predictions_", task_id, ".csv"),
               append = append)

      
      
    # Save into list
    # gam_list[[x]] <- gam_snp_int
    # names(gam_list)[x] <- snp
    
    x = x + 1
    
    
    # Save individual .Rdata file for each model
    # cat("Saving gam model to file... \n")
    # 
    # gam_name <- paste0("gam_", climate_var_dif, "_", snp)
    # assign(gam_name, gam_snp_int) # Assign model the name
    # 
    # save(list = gam_name, file = paste0("./gam_mods_data/", gam_name, ".Rdata"))
    
    # Timing of loop
    end_time <- Sys.time()
    cat("This loop took: ")
    print(end_time - start_time)   
    
  }
  
  
  beep(4)
  
  
  ## Save outlier predictions
    # save(outlier_preds, 
    #   file = paste0("./output/outlier_preds_", Sys.Date(), ".Rdata"))
    

  # ## Plot overlay of outlier SNPS
  # ggplot(outlier_preds, aes(x = tmax_sum_dif_unscaled, y = pred,
  #                           group = factor(snp),
  #                           col = factor(type))) +
  #   geom_vline(aes(xintercept = 0), lty = 2, size = .5) +
  #   geom_vline(xintercept = 4.8, size = .5) +
  #   geom_line(alpha = 0.1,
  #             show.legend = FALSE, size = .15) +
  #   geom_smooth(aes(x = tmax_sum_dif_unscaled, 
  #                   y = pred, group = factor(type)), 
  #                   se = FALSE, size = .5) +
  #  
  #   labs(col = "Genotype") +
  #  # xlim(-2.5, 7.5) +
  #   ylab("Relative growth rate") +
  #   xlab("Tmax transfer distance") +
  #   theme_bw(8) + 
  #   scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5)) +
  #   facet_wrap(~type) +
  #   scale_color_manual(values = c("#FF6633", "grey50", "#6699CC")) + 
  #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #                        panel.grid.minor = element_blank(), 
  #         axis.line = element_line(colour = "black"),
  #         legend.position = "none")
  # 
  # 
  # # Save to file
  #   ggsave(filename = paste0("./figs_tables/Figure 2 - outlier responses ", 
  #                            Sys.Date(), ".pdf"),
  #          units = "cm",
  #          height = 4, width = 6,
  #          useDingbats = FALSE )
    
  
  
  
  # Run summary for each model
  gam_list_summary <- lapply(gam_list, summary)
  
  # Save into dataframe 
  gam_mod_summary_df <- data.frame(
    
    # climate var
    climate_var = climate_var_dif,
    
    # Snp name
    snp = names(gam_list),
    
    # Accession used for prediction
    acc_pred = unlist(lapply(gam_list, function(x) x$accession_pred)),
    
    # Sample size
    n = unlist(lapply(gam_list_summary, function(x) x$n)),
    
    # Sample size of genotypes
    n_gen0 = unlist(lapply(gam_list, function(x) x$n_gen0)),
    n_gen1 = unlist(lapply(gam_list, function(x) x$n_gen1)),
    n_gen2 = unlist(lapply(gam_list, function(x) x$n_gen2)),
    
    # Converged?
    converged = unlist(lapply(gam_list, function(x) x$mgcv.conv)),
    
    # Deviance explained
    dev_explained = unlist(lapply(gam_list_summary, function(x) x$dev.expl)),
  
    # Difference in deviance
    dev_dif = unlist(lapply(gam_list, function(x) x$dev_dif)),
    dev_dif_random = unlist(lapply(gam_list, function(x) x$dev_dif_random)),
    
    # Rsquared adjusted
    rsq_adj = unlist(lapply(gam_list_summary, function(x) x$r.sq)),
  
    # Difference in r squared
    rsq_dif = unlist(lapply(gam_list, function(x) x$rsq_dif)),
    rsq_dif_random = unlist(lapply(gam_list, function(x) x$rsq_dif_random)),
    
    # Estimated degrees of freedom
   # edf = unlist(lapply(gam_list, function(x) sum(x$edf))),


    # Predictions of change in height
    height_change_gen_0 = unlist(lapply(gam_list, function(x) x$height_change_gen_0)),
    height_change_gen_1 = unlist(lapply(gam_list, function(x) x$height_change_gen_1)),
    height_change_gen_2 = unlist(lapply(gam_list, function(x) x$height_change_gen_2)),
  
    # Predictions of absolute height
    height_4.8_gen_0 = unlist(lapply(gam_list, function(x) x$height_4.8_gen_0)),
    height_4.8_gen_1 = unlist(lapply(gam_list, function(x) x$height_4.8_gen_1)),
    height_4.8_gen_2 = unlist(lapply(gam_list, function(x) x$height_4.8_gen_2)),
  
    # Predictions of absolute height
    height_0_gen_0 = unlist(lapply(gam_list, function(x) x$height_0_gen_0)),
    height_0_gen_1 = unlist(lapply(gam_list, function(x) x$height_0_gen_1)),
    height_0_gen_2 = unlist(lapply(gam_list, function(x) x$height_0_gen_2)),
    
    # P values for interaction terms
    p_val_int = ifelse(unlist(lapply(gam_list_summary, function(x)
      length(x$s.table[grep(pattern = "^ti", 
                            rownames(x$s.table)), "p-value"]) == 0)), 
      yes = NA, 
      no = unlist(lapply(gam_list_summary, function(x)
        x$s.table[grep(pattern = "^ti",
                       rownames(x$s.table)), "p-value"]))),
    
    # Use grep to look at the last number of the smooth summary table - 
    # Should correspond to the genotype
    # Returns NA if no match
    p_val_gen_0 = ifelse(unlist(lapply(gam_list_summary, function(x)
                                length(x$s.table[grep(pattern = "0$",
                                               rownames(x$s.table)), "p-value"]) == 0)),
                         yes = NA,
                         no = unlist(lapply(gam_list_summary, function(x)
                                                    x$s.table[grep(pattern = "0$",
                                                   rownames(x$s.table)), "p-value"]))),
    p_val_gen_1 = ifelse(unlist(lapply(gam_list_summary, function(x)
                          length(x$s.table[grep(pattern = "1$",
                            rownames(x$s.table)), "p-value"]) == 0)),
                              yes = NA,
                           no = unlist(lapply(gam_list_summary, function(x)
                                  x$s.table[grep(pattern = "1$",
                                  rownames(x$s.table)), "p-value"]))),

    p_val_gen_2 = ifelse(unlist(lapply(gam_list_summary, function(x)
                           length(x$s.table[grep(pattern = "2$",
                            rownames(x$s.table)), "p-value"]) == 0)),
                          yes = NA,
                          no = unlist(lapply(gam_list_summary, function(x)
                          x$s.table[grep(pattern = "2$",
                                         rownames(x$s.table)), "p-value"])))
    
  )
  
  head(gam_mod_summary_df)

