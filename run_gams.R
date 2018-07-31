# Read in command line arguments
	args<-commandArgs(trailingOnly=TRUE)

  print(args)

	task_id <- as.numeric(args[1])
  interval <- as.numeric(args[2])
  climate_var_dif <- as.character(args[3])

# Load libraries
  library(mgcv)
  library(tidyverse)
  library(visreg) # Installed on home directory
  library(patchwork) # Installed on home directory

# Load data
 load("../gam_cluster_2018-07-31.Rdata")


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
  
  ## Convert fixed and random parts of formula to actual formula objetc
  # No idea why we need to call formula twice, but it works
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
    
    # Make sure there's at least 3 levels
    #  if(length(table(dat_snp[, snp])) < 3){
    #    cat("Skipping because not at least 3 levels")
    #    next
    #  }

    # Choose numbers of knots
    # k = length(table(pull(dat_snp, snp)))
    
    # Convert to factor
    # dat_snp[, snp] <- as.factor(pull(dat_snp, snp))
     dat_snp_unscaled[, snp] <- as.factor(pull(dat_snp_unscaled, snp))
    

# Run gams
   	cat("Working on: ", snp, "... number: ", x, " ...\n" )

  # For SNPS as continuous  
   # fixed_effects_int <- paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"cr\") + te(", snp, ", bs=\"cr\", k = ", k, ") +  te(", climate_var_dif,", bs=\"cr\") + ti(", climate_var_dif,", ",snp,", bs=\"cr\", k =", k," )"))

  # For SNPs as factors
   fixed_effects_int <- paste0(paste0("rgr ~ section_block + ", snp, "+ s(height_2014, bs=\"cr\") +  s(", climate_var_dif,", bs=\"cr\") + s(", climate_var_dif,", by = ", snp, ", bs=\"cr\")")) 
    
    ## Add gen pcs and interactions
  # fixed_effects_int <- paste0(fixed_effects_int, paste0(paste0("+ s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC1_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC2_gen, bs=\"cr\") + s(", climate_var_dif,", by = PC3_gen, bs=\"cr\")")))
    
    ## Add MEMs
  #  fixed_effects_int <- paste0(fixed_effects_int,  "+ s(mem1, bs =\"cr\") + s(mem2, bs =\"cr\") + s(mem3, bs =\"cr\") + s(mem4, bs =\"cr\")")


  # Gam with interaction effect
    gam_snp_int = bam(formula = make_formula(fixed_effects_int, random_effects),
                  data = dat_snp_unscaled[!is.na(pull(dat_snp_unscaled, snp)), ],
                  discrete = TRUE,
                  nthreads = 8,
                  method = "fREML", 
                #  method = "ML",
                  family = "gaussian")

   # Check for convergence
   # if(!gam_snp_int$mgcv.conv){
    #  cat("Model not converged... trying again \n")
    #  gam_snp_int = bam(formula = make_formula(fixed_effects_int, random_effects),
    #                    data = dat_snp, 
    #                    discrete = TRUE,
    #                    nthreads = 8,
    #                    method = "fREML", family = "tw", 
    #                    control = list(trace = FALSE,
    #                                   maxit = 500))
    #  
    #  if(!gam_snp_int$mgcv.conv){
    #    cat("Model STILL not converged...\n")
    #  }
   # }
    

 ## Save sample size per genotype
      gen_table <-  table(dat_snp_unscaled[, snp])
      
      gam_snp_int$n_gen0 <- gen_table["0"]
      gam_snp_int$n_gen1 <- gen_table["1"]
      gam_snp_int$n_gen2 <- gen_table["2"]

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
        dplyr::mutate(!!snp := genotype)

	# Change accession number if genetic data is missing so we don't get errors in prediction
	while(!(pred_sub$accession[1] %in%  gam_snp_int$model$accession)){
	  pred_sub$accession <- as.character(as.numeric(pred_sub$accession) + 1)
	}



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
  #  visreg(gam_snp_int, xvar = climate_var_dif, by = snp,
  #         overlay = TRUE, partial = FALSE, rug = FALSE,
  #         breaks = c(unique(pred_sub$genotype_scaled)),
  #         xtrans = function(x) {(x *  scaled_var_sds_gbs_only[climate_var_dif]) + scaled_var_means_gbs_only[climate_var_dif]},
  #         ylab = "5 year height (cm)")
  #  
  #  dev.off()      


    
  # Make PREDICTION
    gam_preds <-  predict(gam_snp_int, 
                          newdata = pred_sub,
                          type = "response",
                          se = TRUE)
    pred_sub$pred <- gam_preds$fit
    pred_sub$se <- gam_preds$se.fit
    
    pred_sub$genotype <- factor(pred_sub$genotype)
    
   # Backtransform climate variable
     dat_snp_trans <- dat_snp %>%
                  dplyr::select(climate_var_dif, snp, rgr) %>%
                     # Back transform X axis
                     mutate(!!climate_var_dif := back_transform(x = get(climate_var_dif),
                                                     var = climate_var_dif,
                                                     means = scaled_var_means_gbs_only,
                                                     sds = scaled_var_sds_gbs_only)) %>%
                    # Filter down to just 0,1,2 genotypes
                    dplyr::filter(get(snp) %in% c(0, 1, 2))

  
 # # Climate transfer functions by genotype     
 #   p1 <-  ggplot(pred_sub) + 
 #     geom_ribbon(aes(tmax_sum_dif_unscaled, ymin = pred - 2 * se, ymax = pred + 2*se, fill = factor(genotype)), 
 #                 alpha = 0.25, show.legend = FALSE)  +
 #     geom_line(data = pred_sub, aes(x = tmax_sum_dif_unscaled, y = pred, col = factor(genotype)), lwd = 1.5) + 
 #     geom_vline(aes(xintercept = 0), lty = 2) +
 #     labs(col = "Genotype") +
 #     ylab("5 yr height (cm)") +
 #     xlab(climate_var_dif) +
 #     ggtitle(snp) +
 #     theme_bw(15)
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
 #   pdf(paste0("./model_plots/", snp,"_gg.pdf"), width = 15, height = 5)  
 #    print(p1 + p2)
 #   dev.off()
   


   # Make prediction for XXXX degree increase in temperature
    for(genotype in c(0,1,2)){
      
      if(!genotype %in% as.numeric(as.character(pred_sub$genotype))){
        gam_snp_int[paste0("height_change_gen_", genotype)] <- NA
        next
      }
      
      pred_sub_temp <- pred_sub[pred_sub$genotype == genotype, ]
      
      height_change <- (pred_sub_temp$pred[pred_sub_temp$tmax_sum_dif_unscaled == 4.88][1] - 
                          pred_sub_temp$pred[pred_sub_temp$tmax_sum_dif_unscaled == 0][1]) / pred_sub_temp$pred[pred_sub_temp$tmax_sum_dif_unscaled == 4.88][1] * 100 
      
      gam_snp_int[paste0("height_change_gen_", genotype)] <- height_change
      
      }


    # Attach data so that we can use visreg later if we need
    # gam_snp$data <- dat_snp
    

    # Save into list
    gam_list[[x]] <- gam_snp_int
    names(gam_list)[x] <- snp
    
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


  
  
  # Run summary for each model
  gam_list_summary <- lapply(gam_list, summary)
  
  # Save into dataframe 
  gam_mod_summary_df <- data.frame(
    
    # climate var
    climate_var = climate_var_dif,
    
    # Snp name
    snp = names(gam_list),
    
    # Sample size
    n = unlist(lapply(gam_list_summary, function(x) x$n)),

    # Sample size of genotypes
    n_gen0 = unlist(lapply(gam_list, function(x) x$n_gen0)),
    n_gen1 = unlist(lapply(gam_list, function(x) x$n_gen1)),
    n_gen2 = unlist(lapply(gam_list, function(x) x$n_gen2)),

    # Converged?
   # converged = unlist(lapply(gam_list, function(x) x$mgcv.conv)),
    
    # Deviance explained
    dev_explained = unlist(lapply(gam_list_summary, function(x) x$dev.expl)),
    
    # Rsquared adjusted
    r_sq_adj = unlist(lapply(gam_list_summary, function(x) x$r.sq)),
    
    # Estimated degrees of freedom
  #  edf = unlist(lapply(gam_list, function(x) sum(x$edf))),
    
    # AIC
   # aic = unlist(lapply(gam_list, function(x) x$aic)),
    
    # Predictions of height
    height_change_gen_0 = unlist(lapply(gam_list, function(x) x$height_change_gen_0)),
    height_change_gen_1 = unlist(lapply(gam_list, function(x) x$height_change_gen_1)),
    height_change_gen_2 = unlist(lapply(gam_list, function(x) x$height_change_gen_2)),

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
  
 
 write.csv(gam_mod_summary_df, file = paste0("./model_summaries/gam_summaries_", task_id, ".csv"),
                               row.names = FALSE)





