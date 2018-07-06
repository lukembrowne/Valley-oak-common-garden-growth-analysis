

# Load libraries ----------------------------------------------------------
  
  source("./01_clean-process-explore-data.R")
  source("./03_adding-climate-data.R")
  
  dim(dat_all_scaled)
  dim(dat_gbs_only_scaled)


  # install_github("jdstorey/qvalue")
  library(gamm4)
  library(visreg)
  library(rrBLUP)
  library(qvalue)
  library(MuMIn)
  library(vegan)
  library(beepr)
  library(patchwork)
  library(rasterVis)
  
  # function to transform from raw to scaled variables
    forward_transform <- function(x, var, means, sds){
      ( x - means[var]) / sds[var]
    }
    
    # Function to back transform variables
    back_transform <- function(x, var, means, sds){
      (x * sds[var]) + means[var]
    }


# Goal is to run an individual gam for each SNP, like a fancy GWAS, and then correct for multiple testing after

# Format and save to .Rdata for cluster -----------------------------------

  # Joining climate and genetic data
  # Do an inner join to avoid having GBS moms with climate data but poor genotyping data
 
  # Just moms in common garden  
   gen_dat_clim <- inner_join(climate_gbs_mom, 
                             dplyr::select(gen_dat, -accession), 
                             by = c("id" = "gbs_name"))
  dim(gen_dat_clim)
  
  # All GBS samples
  gen_dat_clim_all <- inner_join(climate_gbs_mom, 
                             dplyr::select(gen_dat_all, -accession), 
                             by = c("id" = "gbs_name"))
  dim(gen_dat_clim_all)
  
  
  # Getting genetic data
  snp_col_names <- colnames(gen_dat[, -c(1, 2)])
  
  # Calculate percent missing data  
  missing <- gen_dat_clim %>%
    dplyr::select(id, accession) %>%
    mutate(missing = apply(gen_dat_clim[, snp_col_names], 1, function(x) sum(is.na(x))/length(x)))
  
  missing %>%
    arrange(desc(missing))

  # Scale data and impute missing before converting 012 to -101
  bglr_gen_scaled <- flashpcaR::scale2(gen_dat_clim[, snp_col_names], impute = 2)
  bglr_gen_scaled[1:10, 1:10]
  summary(bglr_gen_scaled[, 1:10])
  
  
  ## Calculate PCA on genetic data
    pca_gen = prcomp(bglr_gen_scaled, center = FALSE, scale = FALSE)
   # biplot(pca_gen)
    summary(pca_gen)
  #  screeplot(pca_gen, bstick = TRUE)
    
    pcs <- as.data.frame(pca_gen$x[, 1:10])
    colnames(pcs) <- paste0(colnames(pcs), "_gen")
    
    bglr_gen_scaled <- bind_cols(bglr_gen_scaled, pcs)

    
  # Set up dataframe with SNP data, height data, and climate data
    dat_snp = left_join(dat_gbs_only_scaled,
                      bind_cols(gen_dat_clim[, c("accession", snp_col_names)],
                                         pcs),
                      by = "accession")
  
   
  # Set up factors  
    dat_snp$accession <- factor(dat_snp$accession)
    dat_snp$section_block <- factor(dat_snp$section_block)
    dat_snp$section <- factor(dat_snp$section)
    
    
  # Scale PCA axes
    dat_snp$PC1_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC1_gen)) / sd(dat_snp$PC1_gen)
    dat_snp$PC2_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC2_gen)) / sd(dat_snp$PC2_gen)
    dat_snp$PC3_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC3_gen)) / sd(dat_snp$PC3_gen)
    dat_snp$PC4_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC4_gen)) / sd(dat_snp$PC4_gen)
    dat_snp$PC5_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC5_gen)) / sd(dat_snp$PC5_gen)
    
    
  # Save unscaled version  
    dat_snp_unscaled <- dat_snp
 
    
  # Scale genotypes  
    scaled_snps_means <- apply(dat_snp[, snp_col_names], MARGIN = 2,
          function(x) mean(x, na.rm = TRUE))
    
    scaled_snps_sds <- apply(dat_snp[, snp_col_names], MARGIN = 2,
                               function(x) sd(x, na.rm = TRUE))
    
    dat_snp[, snp_col_names] <-   apply(dat_snp[, 
                                             snp_col_names], MARGIN = 2,
                    function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    
    dat_snp[1:10, 1050:1060]
    
    # Mean impute
    # dat_snp[, snp_col_names] <- apply(dat_snp[, 
    #                                                 snp_col_names], MARGIN = 2,
    #            function(x) ifelse(is.na(x), yes = 0, no = x))
    
    dat_snp[1:10, 1050:1060]
    
    summary(dat_snp[, 1050:1060])
    
    
    
# Set a prediction data frame


 # Prediction frame with a XXXX degree increase
    
    
    
  degree_increase = 4.88  
    
  # Set up each climate variable individually
  pred <-  expand.grid(height_2014 = 0,
                      accession = "1",
                      section_block = "IFG_1",
                      PC1_gen = 0, PC2_gen = 0, PC3_gen = 0)
  
  for(var in climate_vars_dif){
    
    var_unscaled <- paste0(var, "_unscaled")
    
    # Set scaled climate var dif, with a 0 and 3 degree increase as well
    var_df <- expand.grid(accession = "1",
                          var_temp = c(seq(min(dat_all_scaled[ ,var]),
                                               max(dat_all_scaled[ ,var]),
                                               length.out = 50), 
                                           forward_transform(c(0, degree_increase),
                                                             var = var,
                                                             means = scaled_var_means_gbs_only,
                                                             sds = scaled_var_sds_gbs_only)))
    # Rename columns and set unscaled climate_var dif
    var_df <- var_df %>%
          rename(!!var := var_temp) %>%
          mutate(!!var_unscaled := back_transform(x = get(var),
                                         var = var,
                                         means = scaled_var_means_gbs_only,
                                         sds = scaled_var_sds_gbs_only))
    
   # Do a left join if the first variable, or else just bind cols, removing accession column  
    if(nrow(pred) == 1){
    pred <- left_join(pred, var_df)
    } else {
      pred <- bind_cols(pred, var_df[, -1])
    }
  }
  
  head(pred)

  pred

# Save data to file that will be uploaded to cluster  
  # save(dat_snp, dat_snp_unscaled,
  #      snp_col_names, pred,
  #      scaled_snps_means, scaled_snps_sds,
  #      scaled_var_means_gbs_only, scaled_var_sds_gbs_only,
  #   file = paste0("./output/gam_cluster_", Sys.Date(), ".Rdata"))


  
 
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
  
  
  task_id = 10; interval = 5
  
  
# Loop through snps by index based on task id
  for(snp_index in task_id:(task_id + interval - 1)){
    
    # To avoid going past number of snps
    if(snp_index > length(snp_col_names)){
      next
    }
  
    # Choose snp
    snp <- snp_col_names[snp_index]
  
      # For timing loops
      start_time <- Sys.time()
    
    # Choose numbers of knots
    k = length(table(pull(dat_snp, snp)))
    
    # Make sure there's at least 3 levels
    if(length(table(dat_snp[, snp])) < 3){
      cat("Skipping because not at least 3 levels")
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
    fixed_effects_int <- paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"cr\") + te(", snp, ", bs=\"cr\", k = ", k, ") +  te(", climate_var_dif,", bs=\"cr\") + ti(", climate_var_dif,", ",snp,", bs=\"cr\", k =", k," )"))
    
    fixed_effects_int <- paste0(paste0("rgr ~ section_block + ", snp, "+ s(height_2014, bs=\"cr\") +  s(", climate_var_dif,", bs=\"cr\") + s(", climate_var_dif,", by = ", snp, ", bs=\"cr\")"))
    
   # fixed_effects_int <- paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"ts\") + te(", snp, ", k = ", k, ", bs=\"ts\") +  te(", climate_var_dif,", bs=\"ts\") + ti(", climate_var_dif,", ",snp, ", bs=\"ts\", k = ", k, ")"))
    
    ## Add gen pcs
   # fixed_effects_int <- paste0(fixed_effects_int, "+ s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")")
    
    ## Add MEMs
  #  fixed_effects_int <- paste0(fixed_effects_int,  "+ s(mem1, bs =\"cr\") + s(mem2, bs =\"cr\") + s(mem3, bs =\"cr\") + s(mem4, bs =\"cr\")")
                                       
      # fixed_effects_no_int <- paste0(paste0("height_2017 ~ section_block + s(height_2014, bs=\"cr\") + te(", snp, ", bs=\"cr\", k = ", k, ") +  te(", climate_var_dif,", bs=\"cr\")  + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
    # 
    
    ## Randomize genotype
   # dat_snp[, snp] <- dat_snp[sample(1:nrow(dat_snp)) , snp]
    
    start <- Sys.time()
   # Gam with interaction effect
    gam_snp_int = bam(formula = make_formula(fixed_effects_int, random_effects),
                  data = dat_snp_unscaled, 
               #  discrete = TRUE, # Discrete makes the right curve up!!
                  nthreads = 8,
                method = "fREML", 
              #  method = "ML",
                  family = "gaussian", 
                  control = list(trace = FALSE))
    
    summary(gam_snp_int)

    Sys.time() - start

    ## Look at individuals with outlier residuals
   # resid_outlier <- dat_snp$accession_progeny[residuals(gam_snp_int) < -0.5]
   # View(dat_all[dat_all$accession_progeny %in% resid_outlier, ])

 #  visreg(gam_snp_int, xvar = "tmax_sum_dif", partial = FALSE)

    
    ## Save sample size per genotype
      gen_table <-  table(dat_snp_unscaled[, snp])
      
      gam_snp_int$n_gen0 <- gen_table["0"]
      gam_snp_int$n_gen1 <- gen_table["1"]
      gam_snp_int$n_gen2 <- gen_table["2"]

    # Check for convergence
    # if(!gam_snp_int$mgcv.conv){
    #   cat("Model not converged... trying again \n")
    #   gam_snp_int = bam(formula = make_formula(fixed_effects_int, random_effects),
    #                     data = dat_snp, 
    #                     discrete = TRUE,
    #                     nthreads = 8,
    #                     method = "fREML", family = "tw", 
    #                     control = list(trace = FALSE,
    #                                    maxit = 500))
    #   
    #   if(!gam_snp_int$mgcv.conv){
    #     cat("Model STILL not converged...\n")
    #   }
    #  }

    # gam_snp_no_int =  bam(formula = make_formula(fixed_effects_no_int,
    #                                              random_effects),
    #                       data = dat_snp, 
    #                       discrete = TRUE,
    #                       nthreads = 8,
    #                       method = "fREML", family = "tw")
    
    # Compare AICs
      #gam_snp_int$delta_aic <- gam_snp_int$aic - gam_snp_no_int$aic
      gam_snp_int$delta_aic <- NA
  
    # Visualization plots
    
    # For snps as factors
      # pdf(paste0("./output/model_visualizations/", snp,".pdf"), width = 8, height = 5)
      # visreg(gam_snp_int, xvar = climate_var_dif, 
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
    pred_sub <-  pred_sub  %>%
      dplyr::filter(genotype_scaled %in% pull(dat_snp, snp)) %>%
      dplyr::mutate(!!snp := genotype_scaled)
    
      ## For SNPs as factors
      pred_sub <-  pred_sub  %>%
        dplyr::filter(genotype %in% pull(dat_snp_unscaled, snp)) %>%
        dplyr::mutate(!!snp := genotype)
    
    
  # For SNPs as factors  
    # Change accession number if genetic data is missing so we don't get errors in prediction
    # while(any(is.na(dat_snp[as.character(dat_snp$accession) %in% 
    #                         as.character(pred_sub$accession), snp]))){
    #   pred_sub$accession = as.numeric(pred_sub$accession) + 1
    # }
    
    
    ## Visreg plots
  #  pdf(paste0("./output/model_visualizations/", snp,".pdf"), width = 8, height = 5)

   # With partial residuals
    # visreg(gam_snp_int, xvar = climate_var_dif, by = snp,
    #        overlay = TRUE, partial = TRUE, 
    #        breaks = c(unique(pred_sub$genotype_scaled)),
    #        xtrans = function(x) {(x *  scaled_var_sds_gbs_only[climate_var_dif]) + scaled_var_means_gbs_only[climate_var_dif]},
    #        ylab = "5 year height (cm)")
    
    # # Without partial residuals
    # visreg(gam_snp_int, xvar = climate_var_dif, by = snp,
    #        overlay = TRUE, partial = FALSE, rug = FALSE,
    #        breaks = c(unique(pred_sub$genotype_scaled)),
    #        xtrans = function(x) {(x *  scaled_var_sds_gbs_only[climate_var_dif]) + scaled_var_means_gbs_only[climate_var_dif]},
    #        ylab = "5 year height (cm)")
    
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
    dat_snp_trans <- dat_snp_unscaled %>%
      dplyr::select(climate_var_dif, snp, rgr) %>%
      # Back transform X axis
      mutate(!!climate_var_dif := back_transform(x = get(climate_var_dif),
                                                 var = climate_var_dif,
                                                 means = scaled_var_means_gbs_only,
                                                 sds = scaled_var_sds_gbs_only)) %>%
      # Filter down to just 0,1,2 genotypes
      dplyr::filter(get(snp) %in% c(0, 1, 2))
    
    
    # Climate transfer functions by genotype
      p1 <-  ggplot(pred_sub) +
        geom_ribbon(aes(tmax_sum_dif_unscaled, ymin = pred - 2 * se, ymax = pred + 2*se, fill = factor(genotype)),
                    alpha = 0.25, show.legend = FALSE)  +
        geom_line(data = pred_sub, aes(x = tmax_sum_dif_unscaled, y = pred, col = factor(genotype)), lwd = 1.5) +
        geom_vline(aes(xintercept = 0), lty = 2) +
        geom_vline(xintercept = 4.88, lwd = 1.5) +
        labs(col = "Genotype") +
        ylab("5 yr height (cm)") +
        xlab(climate_var_dif) +
        ggtitle(snp) +
        theme_bw(15)

    # Density plots of genotypes
      p2 <- ggplot(dat_snp_trans,
                   aes(x = get(climate_var_dif), fill = factor(get(snp)))) +
        geom_histogram(binwidth = .5, col = "grey10") +
        geom_vline(aes(xintercept = 0), lty = 2) +
        labs(fill = "Genotype") +
        xlab(climate_var_dif) +
        ggtitle(snp) +
        theme_bw(15)

    # Plot to PDF
      pdf(paste0("./output/model_visualizations/", snp,"_gg.pdf"), width = 15, height = 5)
       print(p1 + p2)
      dev.off()

    
  # Make prediction for 3 degree increase in temperature
    # for(genotype in c(0,1,2)){
    #   
    #   if(!genotype %in% as.numeric(as.character(pred_sub$genotype))){
    #     gam_snp_int[paste0("height_change_gen_", genotype)] <- NA
    #     next
    #   }
    #   
    #   pred_sub_temp <- pred_sub[pred_sub$genotype == genotype, ]
    #   
    #   height_change <- (pred_sub_temp$pred[pred_sub_temp[, climate_var_dif_unscaled] == 3][1] - 
    #                       pred_sub_temp$pred[pred_sub_temp[, climate_var_dif_unscaled] == 0][1]) / pred_sub_temp$pred[pred_sub_temp[, climate_var_dif_unscaled] == 3][1] * 100 
    #   
    #   gam_snp_int[paste0("height_change_gen_", genotype)] <- height_change
    #   
    #   }

    # Attach data so that we can use visreg later if we need
    # gam_snp$data <- dat_snp
    
    # Save into list
    gam_list[[x]] <- gam_snp_int
    names(gam_list)[x] <- snp
    
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
  #  converged = unlist(lapply(gam_list, function(x) x$mgcv.conv)),
    
    # Deviance explained
    dev_explained = unlist(lapply(gam_list_summary, function(x) x$dev.expl)),
    
    # Rsquared adjusted
    r_sq_adj = unlist(lapply(gam_list_summary, function(x) x$r.sq)),
    
    # Estimated degrees of freedom
   # edf = unlist(lapply(gam_list, function(x) sum(x$edf))),
    
    # AIC
   #  aic = unlist(lapply(gam_list, function(x) x$aic)),
    
    # Delta AIC with model without interaction
   # delta_aic = unlist(lapply(gam_list, function(x) x$delta_aic)),
    
    # Predictions of height
    # height_change_gen_0 = unlist(lapply(gam_list, function(x) x$height_change_gen_0)),
    # height_change_gen_1 = unlist(lapply(gam_list, function(x) x$height_change_gen_1)),
    # height_change_gen_2 = unlist(lapply(gam_list, function(x) x$height_change_gen_2)),
    
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
  
  
  
  

# Read in cluster run output ----------------------------------------------

 # path_to_summaries <- "./output/model_summaries_run_3027735_pc_gen/"  # Controlling for structure
  path_to_summaries <- "./output/model_summaries_run_3118140_tmax_sum_dif_rgr_ML/"
 # path_to_summaries <- "./output/model_summaries_run_3130265_tmax_sum_dif_rgr_fREML_discrete/"
 # path_to_summaries <- "./output/model_summaries_run_3131613_tmax_sum_dif_rgr_fREML_not_discrete/"
#  path_to_summaries <- "./output/model_summaries_run_3134990_tmax_sum_dif_rgr_fREML_discrete_random_genotype/"
  
  
  sum_df_raw <- plyr::ldply(list.files(path_to_summaries, full = TRUE), read_csv)
  
  dim(sum_df_raw)
  sum_df_raw
  
  summary(sum_df_raw)
  
  
  ## Missing snps
  snp_col_names[which(!snp_col_names %in% sum_df_raw$snp)]
  
  
  
  

# Sort and arrange data ---------------------------------------------------

  # How many converged?
  # table(sum_df$converged)
  # 
  # # Filter to just those that converged
  # sum_df <- sum_df %>%
  #   dplyr::filter(converged)
  # dim(sum_df)
  # 
  
  ## Filter down to those with just a certain amount of gene copies
  
  gen_copies_thresh = 25
  
  sum_df <- sum_df_raw %>%
    filter(n_gen0 > gen_copies_thresh) %>%
    filter(n_gen1 > gen_copies_thresh) %>%
    filter(n_gen2 > gen_copies_thresh)
  
  dim(sum_df)
  
  
  # Figuring out which allele is beneficial
  sum_df <- sum_df %>%
    mutate(beneficial_gen_0 = ifelse(height_change_gen_0 > 0, 1, 0),
           beneficial_gen_1 = ifelse(height_change_gen_1 > 0, 1, 0),
           beneficial_gen_2 = ifelse(height_change_gen_2 > 0, 1, 0)) %>%
    mutate(beneficial_any = (beneficial_gen_0 == 1 | 
                                  beneficial_gen_1 == 1 |
                                  beneficial_gen_2 == 1))
  
  table(sum_df$beneficial_any)
  table(sum_df$beneficial_gen_0)
  table(sum_df$beneficial_gen_1)
  table(sum_df$beneficial_gen_2)  

  ## Calculating q values  
      # Vignette - https://bioconductor.org/packages/release/bioc/vignettes/qvalue/inst/doc/qvalue.pdf
      
      qvals <- qvalue(p = sum_df %>%
                   #     filter(beneficial_any == 1) %>%
                        dplyr::pull(p_val_int),
                      fdr.level = 0.05)
      
      summary(qvals)
      
      table(qvals$significant)
      
      qvals$lfdr
      
      hist(qvals)
      
      # Assign back into DF
      sum_df <- sum_df %>%
        mutate(q_val = qvals$qvalues)
      
        qvals$pi0 # overall proportion of true null hypotheses
      
      # Sort by lowest q value  
      sum_df %>%
        dplyr::arrange(q_val)
      
      
  ## Using p adjust
    sum_df <-   sum_df %>%
      mutate(p_val_int_adj = p.adjust(p_val_int, method = "fdr"))
    
    summary(sum_df$p_val_int_adj)
    
    table(sum_df$p_val_int_adj <= 0.05)
    table(sum_df$p_val_int_adj <= 0.10)
    
    hist(sum_df$p_val_int_adj, breaks = 20)
    


# Select top snps ---------------------------------------------------------

    # Based on q value
    top_snps = sum_df %>%
      dplyr::filter(q_val < 0.05) %>%
      dplyr::filter(beneficial_gen_0 == 1 | beneficial_gen_1 == 1 | beneficial_gen_2 == 1) %>%
     # dplyr::top_n(-round(nrow(sum_df)*.01), p_val_int) %>% # Take only XXX% of snps
      dplyr::arrange(p_val_int)
    
    top_snps
    
    dim(top_snps)
    
    summary(top_snps)
    
    summary(top_snps$p_val_int)
    
    ## Readd top_snps designation to dataframe
    sum_df$top_snp <- FALSE
    sum_df$top_snp[sum_df$snp %in% top_snps$snp] <- TRUE
    
    
  ## Convert top SNPs genotypes to precense absence
    
    gen_dat_bene <- gen_dat_clim_all # Using all GBS samples
    
    for(snp in top_snps$snp){
      
      sub <- sum_df[sum_df$snp == snp, ]
      
      ben0 <- sub$beneficial_gen_0 == 1 # Check if genotype is beneficial or not
      ben1 <- sub$beneficial_gen_1 == 1
      ben2 <- sub$beneficial_gen_2 == 1
      
     gen_dat_bene[, snp] <-  apply(gen_dat_clim_all[, snp], 1, function(x) {
        if(is.na(x)){ NA
        } else if(x == 0 & is.na(ben0)){ NA
        } else if(x == 0 & ben0){ 1
        } else if(x == 1 & is.na(ben1)){ NA
        } else if(x == 1 & ben1){ 1
        } else if(x == 2 & is.na(ben2)){ NA
        } else if(x == 2 & ben2){ 1
        } else {  0
        }
      }) # End apply
     
    } # End SNP loop
    
    summary(gen_dat_bene[, top_snps$snp])
    
    
  ## For rest of SNPs, convert non-reference alleles to 1, reference alleles to 0 for comparison
    for(snp in snp_col_names[!snp_col_names %in% top_snps$snp]){

      gen_dat_bene[, snp] <-  apply(gen_dat_clim_all[, snp], 1, function(x) {
        if(is.na(x)){ NA
        } else if(x == 0){ 0
        } else if(x == 1){ 1
        } else if(x == 2){ 1
        } else {  0
        }
      }) # End apply
      
    } # End SNP loop
    
    dim(gen_dat_bene)
    
    summary(gen_dat_bene[, snp_col_names[!snp_col_names %in% top_snps$snp][1:50]])
    
    
    
    
  
  ## SCP code to copy over model plots
    for(snp in paste0("snp_", top_snps$snp)[1:15]){
    #  system(paste0("scp lukembro@dtn2.hoffman2.idre.ucla.edu:/u/flashscratch/l/lukembro/qlobata_growth/run_291_tmax_sum_dif/model_plots/", snp, "_gg.pdf ./output/model_visualizations/"))
      
     # gg <-  ggplot(dat_snp, aes(tmax_sum_dif, pull(dat_snp, snp))) +
     #    geom_jitter() +
     #    ylab(snp) +
     #    theme_bw() 
     # 
     # print(gg)
     #  
    #  ggsave(paste0("./output/model_visualizations/",snp, "_scatter.pdf"))
      
    }

# Manhattan plots ---------------------------------------------------------

  library(CMplot)
  
  man_df1 <- data.frame(SNP = snp_col_names, 
                       snp_pos)
  
  man_df <- left_join(man_df1, sum_df[, c("snp", "p_val_int", "top_snp")], by = c("SNP" = "snp"))
  
  man_df$index <- 1:nrow(man_df) ## For x axis plotting
 
  
  head(man_df)
  
  ## Assign color for chromosome
  man_df$color <- man_df$chrom %% 3
  man_df$color[man_df$color == 0] <- "#a6cee3"
  man_df$color[man_df$color == 1] <- "#1f78b4"
  man_df$color[man_df$color == 2] <- "#b2df8a"

  ## Make beneficial alleles bigger and different color
  man_df$cex <- 1
  man_df$cex[man_df$top_snp == TRUE] <- 1.5
  
  man_df$color[man_df$top_snp == TRUE] <- "#33a02c"
  man_df$pch[man_df$top_snp] <- 21 # Outline sig ones
  
  
  
## Custom plot
  
  plot(x = man_df$index, y = -log10(man_df$p_val_int), 
       pch =  20, cex = man_df$cex,
       xlab = "SNP", ylab = "-log10(Pvalue)", 
       las = 1,
       col = man_df$color, # To alternate colors for pseudo-chromosomes
       ylim = c(0, 6),
       bty = "l",
       xaxt = "n")
  
  ## Overplot significant points
  points(x = man_df$index[man_df$top_snp], y = -log10(man_df$p_val_int[man_df$top_snp]),
         pch = 21, col = "black",
         cex = man_df$cex[man_df$top_snp], bg = man_df$color[man_df$top_snp])
  
  ## Choose illustrative non-significant point
  points(x = man_df$index[man_df$chrom == 6 & man_df$p_val_int > 0.8][2],
         y = -log10(man_df$p_val_int[man_df$chrom == 6 & man_df$p_val_int > 0.8][2]),
                    pch = 22, col = "black",
                    cex = 1.5, bg = "black")
  
  ## Add labels
  text(x =  c(man_df$index[man_df$top_snp], 
              man_df$index[man_df$chrom == 6 & man_df$p_val_int > 0.8][2]), 
       y = c(-log10(man_df$p_val_int[man_df$top_snp]), 
             -log10(man_df$p_val_int[man_df$chrom == 6 & man_df$p_val_int > 0.8][2])),
       labels = c(man_df$SNP[man_df$top_snp], 
                  man_df$SNP[man_df$chrom == 6 & man_df$p_val_int > 0.8][2]),
    #   adj = c(1, 0),
       pos = 3,
       cex = 0.6)
  
# Save as pdf
  dev.copy(pdf, paste0("./figs_tables/Figure 4 - manhattan plot interaction term ", 
                       Sys.Date(), ".pdf"),
           width = 8, height = 5)
  dev.off()
  

  
  
  
  
  
  
  # SNP density plot
  # CMplot(man_df, plot.type = "d",
  #        col=c("darkgreen", "yellow", "red"),
  #        file.output = FALSE)
  # 
## Manhattan plot with SNP density  
  dev.off()
  CMplot(man_df, plot.type=c("m"), LOG10=TRUE, threshold=1e-3,
         bin.size=1e6,
         chr.den.col=c("darkgreen", "yellow", "red"), file.output = FALSE)
  
## QQplot  
  dev.off()
  CMplot(man_df, plot.type="q", box=TRUE,file.output = FALSE)
  
# With Q values  
   man_df <- left_join(man_df1, sum_df[, c("snp", "q_val")], by = c("SNP" = "snp"))
   dev.off()
   CMplot(man_df, plot.type=c("m"), LOG10=TRUE, threshold= 0.05,
          bin.size=1e6, ylim = c(0, 2),
          chr.den.col=c("darkgreen", "yellow", "red"), file.output = FALSE)
  

  
   
# Random forest on just beneficial alleles ------------------------
   
   library(randomForest)
   library(AUC)
   library(Boruta) # For variable importance?
   

   ## Convert to factor 
   gf_factor <- gen_dat_bene[ , snp_col_names]
   dim(gf_factor)
   
   ## Count number of positives - where allele is found
   pos_thresh <- 25
   
   # For presences
   ind1 <- which(apply(gf_factor, 2, function(x) sum(x == 1, na.rm = TRUE))
                 < pos_thresh)
   if(length(ind1) > 0) gf_factor <- gf_factor[, -ind1] # Need to check if indexes were found
   
   # For absences
   ind2 <- which(apply(gf_factor, 2, function(x) sum(x == 0, na.rm = TRUE)) < pos_thresh)
   if(length(ind2 ) > 0) gf_factor <- gf_factor[, -ind2]
   
   # # Convert to factor
   # gf_factor <- gf_factor %>%
   #   mutate_all(as.factor) 
   # 
   dim(gf_factor)
   
   
   # TO RANDOMIZE GENOTYPES
   # for(col in 1:ncol(gf_factor)){
   #   gf_factor[, col] <- as.factor(sample(0:1, nrow(gf_factor), replace = TRUE))
   # }
   # summary(gf_factor)
   
 
   
   ## Set up rasters
     library(raster)
   
   
   # Shape outlines
     cali_outline <- readShapePoly("./data/gis/california_outline/california_outline.shp",
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))
     
     lobata_range <- readShapePoly("./data/gis/valley_oak_range/qlobata_refined_grouped.shp",
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
     
     lobata_range_rough <- readShapePoly("./data/gis/valley_oak_range/querloba.shp",
                                         proj4string = CRS("+proj=longlat +datum=WGS84"))
     
    # Elevation map
     dem <- raster("./data/gis/dem/ClimateNA_DEM_cropped.tif")
     
     # Climate rasters
     tmax_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmx1951_1980jja_ave_HST_1513103038/tmx1951_1980jja_ave_HST_1513103038.tif")
     
     tmin_winter <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmn1951_1980djf_ave_HST_1513102910/tmn1951_1980djf_ave_HST_1513102910.tif")
     
     aet <- raster("./data/gis/climate_data/BCM/historical/1951-1980/aet1951_1980_ave_HST_1513103180/aet1951_1980_ave_HST_1513103180.tif")
     bioclim_04 <- raster("./data/gis/climate_data/BCM/historical/1951-1980/bioclim_04.tif")
     bioclim_15 <- raster("./data/gis/climate_data/BCM/historical/1951-1980/bioclim_15.tif")
     bioclim_18 <- raster("./data/gis/climate_data/BCM/historical/1951-1980/bioclim_18.tif")
     bioclim_19 <- raster("./data/gis/climate_data/BCM/historical/1951-1980/bioclim_19.tif")
     
     
    # Function to downscale and project raster 
     process_raster <- function(rast){
       rast <- aggregate(rast, fact = 20) # Make smaller by averaging across cells
       rast_latlon <- projectRaster(rast, crs = CRS("+proj=longlat +datum=WGS84"))
       plot(rast_latlon)
       return(rast_latlon)
     }
     
   # Aggregate and project tmax_raster
     tmax_rast   <- process_raster(tmax_rast)
     tmin_winter <- process_raster(tmin_winter)
     aet         <- process_raster(aet)
     bioclim_04  <- process_raster(bioclim_04)
     bioclim_15  <- process_raster(bioclim_15)
     bioclim_18  <- process_raster(bioclim_18)
     bioclim_19  <- process_raster(bioclim_19)
     
   # Match up resolution and extent of DEM to climate rasters
     dem <- resample(dem, tmax_rast) 
     
     compareRaster(tmax_rast, dem) # Make sure rasters are equal
   
   # Find lat long points
     latlon <- xyFromCell(tmax_rast, 1:ncell(tmax_rast))
     
     rf_rast_df <- data.frame(cell_id = 1:ncell(tmax_rast),
                              longitude = latlon[, 1],
                              latitude = latlon[, 2],
                              tmax_sum = raster::extract(tmax_rast, 1:ncell(tmax_rast)),
                              tmin_winter = raster::extract(tmin_winter, 1:ncell(tmin_winter)),
                              aet = raster::extract(aet, 1:ncell(aet)),
                              elevation = raster::extract(dem, 1:ncell(dem)),
                              bioclim_04 = raster::extract(bioclim_04, 1:ncell(bioclim_04)),
                              bioclim_15 = raster::extract(bioclim_15, 1:ncell(bioclim_15)),
                              bioclim_18 = raster::extract(bioclim_18, 1:ncell(bioclim_18)),
                              bioclim_19 = raster::extract(bioclim_19, 1:ncell(bioclim_19)),
                              random = rnorm(1:ncell(dem)))
     rf_rast_df
     
     dim(rf_rast_df)
     
     # Remove NAs from dataframe or random forest will complain
     rf_rast_df_nona <-  rf_rast_df %>%
       dplyr::filter(!is.na(tmax_sum)) %>%
       dplyr::filter(!is.na(elevation)) %>%
       dplyr::filter(!is.na(aet)) %>%
       dplyr::filter(!is.na(bioclim_04))
     
     summary(rf_rast_df_nona)
     dim(rf_rast_df_nona)
   
     
     
   # TO COMPLETELY RANDOMIZE ENVIRONMENTAL VARIABLES
   
     
   
   # for(var in climate_vars_rf){
   #   gen_dat_ben_randomized[, var] <- rnorm(length(pull(gen_dat_ben_randomized, var)),
   #                                          mean = mean(pull(gen_dat_ben_randomized, var)),
   #                                          sd = sd(pull(gen_dat_ben_randomized, var)))
   # }
   # 
     
  # Choose climate variables   
     climate_vars_rf <-  c("tmax_sum", 
                           "tmin_winter", 
                           # "cwd",
                           "aet",
                           "bioclim_04",
                           "bioclim_15",
                           "bioclim_18",
                           "bioclim_19",
                           "elevation",
                           "latitude",
                           "longitude",
                           "random")
     
     pairs.panels(gen_dat_bene[, climate_vars_rf])
     
     
# Initialize variables
     
     stack_top <- stack() # For Top SNPs
     stack_random_snps <- stack() # For Random SNPs
     stack_random_env <- stack() # For randomizing environmental variables

     start_flag = FALSE
     
     reps = 20
     
     n_snps_to_sample <- 100
 
     
## Begin Random forest loop         
  for(mode in c("top_snps", "random_snps")){   #random_env
    
    for(rep in 1:reps){  
      
      # Don't do reps of random_snps
      if(mode == "random_snps" & rep > 1){
        break
      }
    
    # Choose SNPS to loop over  - same SNPs for either top snps or random_env
      if(mode == "top_snps" | mode == "random_env") {
        snps <- top_snps$snp[top_snps$snp %in% colnames(gf_factor)]
      }
      
      # If random_snps, choose same number of top snps, but snps that aren't top snps
      if(mode == "random_snps") {
        snps <- sample(x = colnames(gf_factor)[!colnames(gf_factor) %in% top_snps$snp],
                                          size = n_snps_to_sample)
      }

      # Loop over SNPs
       for(snp in 1:length(snps)){
         
         snp_name <- snps[snp]
         
         cat("Working on mode:", mode, "rep: ", rep, " snp:", snp, "  ", snp_name, "...\n")
         
         # DO NOT RANDOMIZE CLIMATE VARIABLES
          X = gen_dat_bene[, climate_vars_rf]
         
         ## IF mode == random_env permute CLIMATE VARIABLES 
         if(mode == "random_env") {
           for(col in 1:ncol(X)){
             X[, col] <- X[sample(1:nrow(gen_dat_bene)), col]
           }
           
         }

         y = pull(gf_factor, snp_name)
         n1 <- sum(y == 1, na.rm = TRUE)
         n0 <- sum(y == 0, na.rm = TRUE)
         
         not_NAs <- !is.na(y) # True if not NA
         
         # Convert to factor
         y <- factor(y)
         
         # # Run random forest model
         # rf <- randomForest(y = y[not_NAs], # Don't include NAs
         #                      x = X[not_NAs, ], 
         #                      importance = TRUE,
         #                      nodesize = 10,
         #                    sampsize = c(min(n0, n1), min(n0, n1)),
         #                   # sampsize = c(50, min(n0, n1)),
         #                      #  replace = FALSE,
         #                   #   classwt = c(.5, .5),
         #               #   strata = y[not_NAs],
         #                      ntree = 1000)
         # 
         # 
         # rf
         
        
         
         
         dat_rf <- cbind(y = y[not_NAs], X[not_NAs, ])
         
         
         
         boruta_rf = Boruta(#x =  X[not_NAs, ], y = y[not_NAs],
                       y ~ ., data = dat_rf,
                       min.node.size = 10,
                       ntree = 1000,
                       class.weights = c(min(n0, n1), min(n0, n1)),
                       num.threads = 6)
         boruta_rf
         
         attStats(boruta_rf)
         
         
         # Run RF with only significant predictors
         if(any(attStats(boruta_rf)$decision == "Confirmed")){
         
             rf <- randomForest(getConfirmedFormula(boruta_rf),
                                data = dat_rf,
                                importance = TRUE,
                                nodesize = 10,
                                sampsize = c(min(n0, n1), min(n0, n1)),
                                # sampsize = c(50, min(n0, n1)),
                                #  replace = FALSE,
                                #   classwt = c(.5, .5),
                                #   strata = y[not_NAs],
                                ntree = 1000)
             
             
             rf
         } else {
           rf <- NA
         }
         
         # plot(test)
         # plotImpHistory(test)
         
         # print(importance(rf))
         #  cat("OOB error rate: ", mean(rf$err.rate[, 1]), "\n")
        # cat("Class1 error rate: ", rf$confusion[2,3], "\n")
         
         
         # ROC plot
         # # The closer the curve follows the left-hand border and then the top border of the ROC space, the more accurate the test. The closer the curve comes to the 45-degree diagonal of the ROC space, the less accurate the test.
         # roc0 <- roc(rf$votes[,1],factor(1 * (rf$y==0)))
         # roc1 <- roc(rf$votes[,2],factor(1 * (rf$y==1)))
         # 
         # plot(roc0,
         #      main= snp_name, col = "blue", lwd = 2)
         # plot(roc1, add = TRUE, col = "green", lwd = 2)
         # legend("topright", legend = c(0, 1), col = c("blue", "green"), lty = 1, lwd = 2)
         # 
         # 
         
         #  plot(rf)
         
         #  varImpPlot(rf)
         
         # Save importance and OOB error rates
         if(start_flag == FALSE){
           
          # # Importance 
          #  importance_0 <- data.frame(mode = mode, rep = rep, snp = snp_name, type = "gen0", t(randomForest::importance(rf)[, 1]))
          #  importance_1 <- data.frame(mode = mode, rep = rep, snp = snp_name, type = "gen1", t(randomForest::importance(rf)[, 2]))
          #  importance_acc <- data.frame(mode = mode, rep = rep, snp = snp_name, type = "acc", t(randomForest::importance(rf)[, 3]))
          #  importance_gini <- data.frame(mode = mode, rep = rep, snp = snp_name, type = "gini", t(randomForest::importance(rf)[, 4]))
          #  
          #  importance_df <- rbind(importance_0, importance_1, importance_acc, importance_gini)
          #  
           importance_boruta <- data.frame(mode = mode, rep = rep, snp = snp_name,
                                           climate_var = rownames(attStats(boruta_rf)),
                                           meanImp = attStats(boruta_rf)$meanImp,
                                           decision = attStats(boruta_rf)$decision)
           
          # Error rates
        #   error_df <- data.frame(mode = mode, rep = rep, snp = snp_name, tail(rf$err.rate, 1)) # Last tree error rates
           
         } else { 
          # importance_df <- rbind(importance_df, data.frame(mode = mode, rep = rep, snp = snp_name, type = "gen0", t(randomForest::importance(rf)[, 1])))
          # importance_df <- rbind(importance_df, data.frame(mode = mode, rep = rep, snp = snp_name, type = "gen1", t(randomForest::importance(rf)[, 2])))
          # importance_df <- rbind(importance_df, data.frame(mode = mode, rep = rep, snp = snp_name, type = "acc", t(randomForest::importance(rf)[, 3])))
          # importance_df <- rbind(importance_df, data.frame(mode = mode, rep = rep, snp = snp_name, type = "gini", t(randomForest::importance(rf)[, 4])))
          # 
           
          importance_boruta <- rbind(importance_boruta, data.frame(mode = mode, rep = rep, snp = snp_name,
                                                                   climate_var = rownames(attStats(boruta_rf)),
                                                                   meanImp = attStats(boruta_rf)$meanImp,
                                                                   decision = attStats(boruta_rf)$decision))
          
        #  error_df <- rbind(error_df, data.frame(mode = mode, rep = rep, snp = snp_name, tail(rf$err.rate, 1)))
         }
  
         # Plot partial plot
         
         # pp_out <- partialPlot(rf, pred.data= X, 
         #                       x.var = "tmax_sum", plot = TRUE)
         # pp_tmax_temp <- data.frame(snp = snp_name,
         #                       tmax_prob = pp_out$y,
         #                       tmax_sum = pp_out$x)
         # 
         # 
         # if(x == 1){
         #   pp_tmax = pp_tmax_temp
         # } else {
         #   pp_tmax <- bind_rows(pp_tmax, pp_tmax_temp)
         # }
         # 
         
         #   
         # }
         
         
         ## Predict across a raster and add to stack
  
             # Second column is positives
         
         # If no RF with confirmed variables set predictions as NA
         if(is.na(rf)){
           
           rf_rast_df_nona$pred <- NA
           
         } else {
         
             rf_rast_df_nona$pred = predict(rf, rf_rast_df_nona, type = "prob")[, 2]
             
         }

             rf_rast_df_temp <- left_join(rf_rast_df, rf_rast_df_nona[, c("cell_id", "pred")])

             rast <- tmax_rast
             values(rast) <- rf_rast_df_temp$pred
             names(rast) <- snp_name # Name the raster

             #  levelplot(rast, margin = FALSE)
          if(mode == "top_snps") stack_top <- stack(stack_top, rast)
          if(mode == "random_snps")   stack_random_snps <- stack(stack_random_snps, rast)
          if(mode == "random_env")   stack_random_env <- stack(stack_random_env, rast)
        
           start_flag <- TRUE # Flip flag to true
  
       } # End loop over SNPs
    
     } # End rep
    
  } # End mode loop
     
   ## Errors
     table(error_df$rep, error_df$mode)
     
     # save(importance_boruta, stack_top, stack_random_snps,
     #      file = "./output/random_snps_stack.Rdata")
     
     
     error_df_gg <- error_df %>%
       rename(gen0_error = X0, gen1_error = X1) %>%
       gather(key = "error_type", value = "error", -mode, -snp, -rep) %>%
       group_by(mode, snp, error_type) %>% # Average across reps by snp
       summarise_all(mean)
     
     ggplot(error_df_gg, aes(error_type, error, fill = mode)) + 
       geom_boxplot() + theme_bw(15) + geom_hline(yintercept = 0.5, lty = 2) 
   
   
   ## Importance
     
     head(importance_boruta)
     
     importance_boruta$decision <- as.character(importance_boruta$decision)
     importance_boruta$decision[importance_boruta$decision == "Confirmed"] <- 1
     importance_boruta$decision[importance_boruta$decision == "Tentative"] <- 0
     importance_boruta$decision[importance_boruta$decision == "Rejected"] <- 0
     importance_boruta$decision <- as.numeric(importance_boruta$decision)
     table(importance_boruta$decision)
     
    ## Boruta importance 
     importance_boruta_df_top <- importance_boruta %>%
       group_by(mode, snp, climate_var ) %>%
       summarise_all(mean) %>%
       dplyr::select(-snp, -rep) %>%
       group_by(mode, climate_var) %>%
       summarise_all(mean) %>%
       arrange(mode, meanImp) %>%
       filter(mode == "top_snps") %>%
       ungroup() %>%
       dplyr::select(-meanImp, -decision, -mode, -snp) %>%
       mutate(order = row_number()) 
     
     
     importance_boruta_means <- 
                 importance_boruta %>%
                 group_by(mode, snp, climate_var) %>%  ## Average across reps by snp
                 summarise_all(mean) %>%
                 dplyr::select(-snp, -rep) %>%
                 group_by(mode, climate_var) %>%
                 summarise_all(mean)  %>%
                 left_join(., importance_boruta_df_top, by = c("climate_var"))
     
     
     # Make ggplot
     ggplot(importance_boruta_means, aes(order, meanImp, fill = mode)) + 
      geom_bar(stat = "identity", col = "black", 
               position = position_dodge2(reverse = TRUE)) + 
       theme_bw(15) + 
       #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
       ylab("Variable importance") + xlab("") +
       scale_x_continuous(
         breaks = importance_means$order,
         labels = importance_means$climate_var, expand = c(0,0)) + 
       coord_flip()
     
     
     ggplot(importance_boruta_means, aes(order, decision, fill = mode)) + 
       geom_bar(stat = "identity", col = "black", position = position_dodge2(reverse = TRUE)) + 
       theme_bw(15) + 
       #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
       ylab("Proportion significant") + xlab("") +
       scale_x_continuous(
         breaks = importance_means$order,
         labels = importance_means$climate_var, expand = c(0,0)) + 
       coord_flip()
     
     
     
     
    # Order bars by top snps 
   importance_means_only_top <-   importance_df %>%
       gather(key = "climate_var", value = "importance", -snp, -type, -mode, -rep) %>%
      group_by(mode, snp, type, climate_var) %>%  ## Average across reps by snp
      summarise_all(mean) %>%
      dplyr::select(-snp, -rep) %>%
      group_by(mode, type, climate_var) %>%
     summarise_all(mean) %>%
     arrange(mode, type, importance) %>%
     filter(mode == "top_snps") %>%
     ungroup() %>%
     dplyr::select(-importance, -mode, -snp) %>%
     mutate(order = row_number()) 
   
   # Combine into big DF with random snps too
   importance_means <- 
     importance_df %>%
     gather(key = "climate_var", value = "importance", -snp, -type, -mode, -rep) %>%
     group_by(mode, snp, type, climate_var) %>%  ## Average across reps by snp
     summarise_all(mean) %>%
     dplyr::select(-snp, -rep) %>%
     group_by(mode, type, climate_var) %>%
     summarise_all(mean)  %>%
     left_join(., importance_means_only_top, by = c("type", "climate_var"))
   
   importance_means
  table(importance_means$order) # Should be multiple in each order
     
  # Make ggplot
     ggplot(importance_means, aes(order, importance, fill = mode)) + 
       geom_bar(stat = "identity", col = "black", position = position_dodge2(reverse = TRUE)) + 
       theme_bw(15) + facet_wrap(~ type, scales = "free", ncol = 2) + 
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
       ylab("Variable importance") + xlab("") +
       scale_x_continuous(
         breaks = importance_means$order,
         labels = importance_means$climate_var, expand = c(0,0)) + 
       coord_flip()
     
   
 # Looking at rasters
    stack_top
    
    stack_random
     
    ## Take mean of rasters, weighted by class1 error - higher weights to lower error
    # top_snps_stack <- weighted.mean(stack_top,
    #                                 w = 1 - error_df$OOB[error_df$mode == "top_snps"])
    # random_snps_stack <- weighted.mean(stack_random_snps,
    #                                    w = 1 - error_df$OOB[error_df$mode == "random_snps"])
    # random_env_stack <- weighted.mean(stack_random_env, 
    #                                   w = 1 - error_df$OOB[error_df$mode == "random_env"]) 
    
    top_snps_stack <- mean(stack_top)
    random_snps_stack <- mean(stack_random_snps, na.rm = TRUE)
    
 
     levelplot(top_snps_stack, margin = FALSE,
               main = "Probability of finding beneficial allele") +
       latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey10")) +
       latticeExtra::layer(sp.polygons(lobata_range_rough, col = "black")) +
       latticeExtra::layer(sp.polygons(lobata_range, col = "black"))
     
   # Get values by region
   # Not 100% sure if extract is doing exactly what we want   
    prob_by_region <- data.frame(region = lobata_range$REGION,
                                 prob   = unlist(lapply(raster::extract(top_snps_stack,
                                                          y = lobata_range), mean)))
    
    ggplot(prob_by_region, aes(region, prob)) + 
      geom_bar(stat = "identity", col = "black", fill = "steelblue2") + theme_bw(15)
     
     
     # Random SNPs
     levelplot(random_snps_stack, margin = FALSE,
               main = "Random Snps")
     
     # Randomize environmental variables
     levelplot(random_env_stack, margin = FALSE,
               main = "Randomdomized environmental variables")
     
     
     # Difference between beneficial and random SNPs
     # Set a 0 point 
     levelplot(top_snps_stack - random_snps_stack, margin = FALSE, 
               main = c("Dif bw bene and random SNP")) +
       latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey10")) +
       latticeExtra::layer(sp.polygons(lobata_range_rough, col = "black"))
     
     
     levelplot(top_snps_stack - random_env_stack, margin = FALSE, 
               main = c("Dif bw bene and random environmental variables")) +
       latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey10")) +
       latticeExtra::layer(sp.polygons(lobata_range_rough, col = "black"))
     
   
     
     
     
     
   
   # Plot partial dependence plots of TMax
   ggplot(pp_tmax, aes(x = tmax_sum, y = tmax_prob,
                       group = snp, alpha = .25)) + geom_line() + theme_bw()
   
   
   
   
 ## Calculate correlation between prob of finding beneficial allele and future climate change
   
   # Read in future tmax scenarios
   tmax_CCSM4 <- raster("./data/gis/climate_data/BCM/future/CCSM4_rcp85/tmx2070_2099jja_ave_CCSM4_rcp85_1529358915/tmx2070_2099jja_ave_CCSM4_rcp85_1529358915.tif")
   tmax_CCSM4 <- projectRaster(tmax_CCSM4, crs = CRS("+proj=longlat +datum=WGS84"))
   tmax_CCSM4 <- resample(tmax_CCSM4, tmax_rast)
   compareRaster(tmax_CCSM4, tmax_rast)
   levelplot(tmax_CCSM4, margin = FALSE)
   
   
   tmax_CNRM <- raster("./data/gis/climate_data/BCM/future/CNRM_rcp85/tmx2070_2099jja_ave_CNRM_rcp85_1529358976/tmx2070_2099jja_ave_CNRM_rcp85_1529358976.tif")
   tmax_CNRM <- projectRaster(tmax_CNRM, crs = CRS("+proj=longlat +datum=WGS84"))
   tmax_CNRM <- resample(tmax_CNRM, tmax_rast)
   compareRaster(tmax_CNRM, tmax_rast)
   levelplot(tmax_CNRM, margin = FALSE)
   
   tmax_IPSL <- raster("./data/gis/climate_data/BCM/future/IPSL_rcp85/tmx2070_2099jja_ave_IPSL_rcp85_1529358959/tmx2070_2099jja_ave_IPSL_rcp85_1529358959.tif")
   tmax_IPSL <- projectRaster(tmax_IPSL, crs = CRS("+proj=longlat +datum=WGS84"))
   tmax_IPSL <- resample(tmax_IPSL, tmax_rast)
   compareRaster(tmax_IPSL, tmax_rast)
   levelplot(tmax_IPSL, margin = FALSE)
   
   tmax_MIROC <- raster("./data/gis/climate_data/BCM/future/MIROC_rcp85/tmx2070_2099jja_ave_MIROC_rcp85_1529357403/tmx2070_2099jja_ave_MIROC_rcp85_1529357403.tif")
   tmax_MIROC <- projectRaster(tmax_MIROC, crs = CRS("+proj=longlat +datum=WGS84"))
   tmax_MIROC <- resample(tmax_MIROC, tmax_rast)
   compareRaster(tmax_MIROC, tmax_rast)
   levelplot(tmax_MIROC, margin = FALSE)
   
   # Stack rasters
    future_stack <- stack(tmax_CCSM4, tmax_CNRM, tmax_IPSL, tmax_MIROC)
  
    levelplot(future_stack, margin = FALSE)
   
   # Calculate future differences in temperature
    future_stack_dif <- future_stack - tmax_rast 

    future_stack_dif <- mean(future_stack_dif)
    levelplot(future_stack_dif, margin = FALSE)
    
    
    # Correlations in future temp differences and probability of finding allele
      data.frame(tmax_sum_dif = values(future_stack_dif),
                 prob_bene = values(bene_stack)) %>%
        ggplot(., aes(tmax_sum_dif, prob_bene)) + geom_point(size = 0.15, alpha = 0.25) +
        geom_smooth(method = "lm") + theme_bw(15) +
        geom_smooth(col = "forestgreen") 
      
    cor(values(future_stack_dif), values(bene_stack),
        use = "na.or.complete")
    
    
    # Just lobata range
      future_stack_dif_range <- mask(future_stack_dif, lobata_range_rough)
      
      plot(future_stack_dif_range, bene_stack_range,
           las = 1, ylab = "Probability of finding beneficial allele",
           xlab = "Tmax_sum_dif")
      
      cor(values(future_stack_dif_range), values(bene_stack_range),
          use = "na.or.complete")
      
    data.frame(tmax_sum_dif = values(future_stack_dif_range),
                         prob_bene = values(bene_stack_range)) %>%
      ggplot(., aes(tmax_sum_dif, prob_bene)) + geom_point(size = 0.15, alpha = 0.25) +
        geom_smooth(method = "lm") + geom_smooth(col = "forestgreen") + theme_bw(15)
  
      
      
      rasterCo  
   
      library(spatialEco)
      
     test <-  rasterCorrelation(future_stack_dif, bene_stack, s = 7)
     test
     levelplot(test, margin = FALSE)
   
   


  

# Gradient forest test ----------------------------------------------------
  
  library(gradientForest)
  
  climate_vars_gf <- c("tmax_sum", 
                      # "tmax_sum_lgm",
                      # "tmin_winter", 
                     # "bioclim_04",
                       #"tmin_winter_lgm", 
                       #  "DD5",
                     #  "DD5_lgm",
                       "random", 
                      "latitude",
                     # "mem1",
                    #  "mem2",
                    #  "mem3",
                    #  "mem4",
                      "longitude", 
                      "elevation")
  
  pairs.panels(gen_dat_clim[, climate_vars_gf])
  
  
  ## RANDOM SUBSAMPLE OF SNPS
  # gf_sample <- sample(snp_col_names, length(top_snps$snp))
  
  
  gf_scaled <- flashpcaR::scale2(gen_dat_clim[, gf_sample], impute = 2)
  
  
  gen_dat_clim_randomized <- gen_dat_clim[, climate_vars_gf]
  

  ## Factors
  gf_factor <- gen_dat_bene[ ,top_snps$snp]
  gf_factor <- gf_factor %>%
                 replace(is.na(.), 0) %>%
                mutate_all(as.factor) 
  
  ## Count number of positives
  pos_thresh <- 20
  gf_factor <- gf_factor[, -which(apply(gf_factor, 2, function(x) sum(x == 1)) < pos_thresh)]
  gf_factor <- gf_factor[, -which(apply(gf_factor, 2, function(x) sum(x == 0)) < pos_thresh)]
  
  dim(gf_factor)
  summary(gf_factor)
              
  
  x = 1
  r2_ran <- NULL
  num_pos <- NULL
  
  
  maxLevel <- log2(0.368*nrow(gen_dat_clim_randomized)/2) #account for correlations, see ?gradientForest 
  
  
  while(x < 10){
    
    cat("working on:", x, "...\n")
    
  #  gen_dat_clim_randomized <- gen_dat_clim_randomized[sample(1:nrow(gen_dat_clim_randomized)),]
    
      gf <- gradientForest(cbind(gen_dat_clim_randomized,
                                       gf_factor ),
                           predictor.vars = climate_vars_gf, 
                          response.vars = colnames(gf_factor),
                         #response.vars = gsub("snp_", "", sum_df$snp),
                          trace = TRUE,
                         ntree = 100, mtry = 3,
                         maxLevel = maxLevel, # From Fitzpatrick R script
                         corr.threshold = 0.5)
  
      gf$result # Rsquared of positive loci
      gf$species.pos.rsq
      
      cat("\n", mean(gf$result), "\n")
      
      importance(gf)
      
      # IF NULL
      if(length(gf) == 0){
        if(x == 1) {next}
        num_pos[x] <- 0
        r2_ran[x] <- 0
        importance_df <- rbind(importance_df, 0)
        x = x+1
        next
      }
      
      num_pos[x] <- gf$species.pos.rsq
      r2_ran[x] <- mean(gf$result)
      if(x == 1) {importance_df <- importance(gf)}
      importance_df <- rbind(importance_df, importance(gf))
      x = x +1
  }
  
  r2_ran
  mean(r2_ran)
  
  num_pos
  mean(num_pos)
  
  importance_df
  sort(colMeans(importance_df), decreasing = TRUE)
  
  plot(gf, plot.type = "O")
  
  
  ## Need to make rasters with each environmental variable
  ## And then extract this data to a dataframe to use for predictions in GF
  
  library(raster)
  
  tmax_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmx1951_1980jja_ave_HST_1513103038/tmx1951_1980jja_ave_HST_1513103038.tif")
  tmin_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmn1951_1980djf_ave_HST_1513102910/tmn1951_1980djf_ave_HST_1513102910.tif")
  
  process_raster <- function(rast){
    rast_temp <- aggregate(rast, fact = 10) # Make smaller
    rast_latlon <- projectRaster(rast_temp, crs = CRS("+proj=longlat +datum=WGS84"))
    plot(rast_latlon)
    return(rast_latlon)
  }
  
  tmax_rast <- process_raster(tmax_rast)
  tmin_rast <- process_raster(tmin_rast)
  
  # Find lat long points
  latlon <- xyFromCell(tmax_rast, 1:ncell(tmax_rast))
  
  gf_rast_df <- data.frame(cell_id = 1:ncell(tmax_rast),
                        longitude = latlon[, 1],
                        latitude = latlon[, 2],
                        tmax_sum = raster::extract(tmax_rast, 1:ncell(tmax_rast)),
                        tmin_winter = raster::extract(tmin_rast, 1:ncell(tmin_rast)))
  gf_rast_df
  
  gf_rast_nona <-  gf_rast_df %>%
    filter(!is.na(tmax_sum))

  
  # transform env using gf models, see ?predict.gradientForest
  gf_pred <- predict(gf, gf_rast_nona[, colnames(gf_rast_nona) %in% climate_vars_gf]) # remove cell column before transforming
  
  
  # map continuous variation - reference SNPs
  refRGBmap <- pcaToRaster(gf_pred, gf_rast_df, tmax_rast, gf_rast_nona$cell_id)
  plotRGB(refRGBmap)

  # Mapping spatial genetic variation --------------------------------------------
  ###### functions to support mapping #####
  # builds RGB raster from transformed environment
  # snpPreds = dataframe of transformed variables from gf or gdm model
  # rast = a raster mask to which RGB values are to be mapped
  # cellNums = cell IDs to which RGB values should be assigned
  pcaToRaster <- function(gf_pred, gf_pred_wna, rast, mapCells){
    require(raster)
    
    pca <- prcomp(gf_pred, center=TRUE, scale.=FALSE)
    
    ##assigns to colors, edit as needed to maximize color contrast, etc.
    a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
    r <- a1+a2; g <- -a2; b <- a3+a2-a1
    
    ##scales colors
    scalR <- (r-min(r))/(max(r)-min(r))*255
    scalG <- (g-min(g))/(max(g)-min(g))*255
    scalB <- (b-min(b))/(max(b)-min(b))*255
    
    ## Join back into data frame that includes NAs
    temp <- data.frame(cell_id = mapCells,
               scalR = scalR,
               scalG = scalG, 
               scalB = scalB)
    
    temp_wna <- left_join(gf_pred_wna, temp, by = "cell_id")
    
    
    ## Check to make sure lengths match up
    if(nrow(temp_wna) != ncell(rast)){
      stop("Error: length of raster and of PCA predictions do not match")
    }
    
    ##assigns color to raster
    rast1 <- rast2 <- rast3 <- rast
    values(rast1) <- temp_wna$scalR
    values(rast2) <- temp_wna$scalG
    values(rast3) <- temp_wna$scalB
    ##stacks color rasters
    outRast <- stack(rast1, rast2, rast3)
    return(outRast)
  }
  
  # Function to map difference between spatial genetic predictions
  # predMap1 = dataframe of transformed variables from gf or gdm model for first set of SNPs
  # predMap2 = dataframe of transformed variables from gf or gdm model for second set of SNPs
  # rast = a raster mask to which Procrustes residuals are to be mapped
  # mapCells = cell IDs to which Procrustes residuals values should be assigned
  RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
    require(vegan)
    PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
    PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
    diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
    residMap <- residuals(diffProcrust)
    rast[mapCells] <- residMap
    return(list(max(residMap), rast))
  }
  
  # OK, on to mapping. Script assumes:
  # (1) a dataframe named env_trns containing extracted raster data (w/ cell IDs)
  # and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
  #
  # (2) a raster mask of the study region to which the RGB data will be written
  
  # transform env using gf models, see ?predict.gradientForest
  predRef <- predict(gfRef, env_trns[,-1]) # remove cell column before transforming
  predGI5 <- predict(gfGI5, env_trns[,-1])
  
  # map continuous variation - reference SNPs
  refRGBmap <- pcaToRaster(predRef, mask, env_trns$cell)
  plotRGB(refRGBmap)
  writeRaster(refRGBmap, "/.../refSNPs_map.tif", format="GTiff", overwrite=TRUE)
  
  # map continuous variation - GI5 SNPs
  GI5RGBmap <- pcaToRaster(predGI5, mask, env_trns$cell)
  plotRGB(refRGBmap)
  writeRaster(refRGBmap, "/.../GI5SNPs_map.tif", format="GTiff", overwrite=TRUE)
  
  # Difference between maps (GI5 and reference) 
  diffGI5 <- RGBdiffMap(predRef, predGI5, rast=mask, mapCells=env_trns$cell)
  plot(diffGI5[[2]])
  writeRaster(diffGI5[[2]], "/.../diffRef_GI5.tif", format="GTiff", overwrite=TRUE)
  ################################################################################
  
  
  # Calculate and map "genetic offset" under climate change ----------------------
  # Script assumes:
  # (1) a dataframe of transformed env. variables for CURRENT climate 
  # (e.g., predGI5 from above).
  #
  # (2) a dataframe named env_trns_future containing extracted raster data of 
  # env. variables for FUTURE a climate scenario, same structure as env_trns
  
  # first transform FUTURE env. variables
  projGI5 <- predict(gfGI5, env_trns_future[,-1])
  
  # calculate euclidean distance between current and future genetic spaces  
  genOffsetGI5 <- sqrt((projGI5[,1]-predGI5[,1])^2+(projGI5[,2]-predGI5[,2])^2
                       +(projGI5[,3]-predGI5[,3])^2+(projGI5[,4]-predGI5[,4])^2
                       +(projGI5[,5]-predGI5[,5])^2+(projGI5[,6]-predGI5[,6])^2
                       +(projGI5[,7]-predGI5[,7])^2)
  
  # assign values to raster - can be tricky if current/future climate
  # rasters are not identical in terms of # cells, extent, etc.
  mask[env_trns_future$cell] <- genOffsetGI5
  plot(mask)
  
  
  
  ################################################################################
  # END SCRIPTS FOR GRADIENT FOREST MODELING
  ################################################################################
  
  
  
  
  
  
  
  
  
  
    plot(gf, plot.type = "O")
    
    most_important <- names(importance(gf))
    
    pred <- expand.grid(latitude = seq(31, 45, by = .1))
    pred$pred <- predict(gf, newdata = pred)
    
    pred
    
    plot(pred$latitude, pred$pred$latitude)
  
    
    plot(gf, plot.type = "S", imp.vars = most_important,
          leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
          cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
          0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
    
    plot(gf, plot.type = "C", imp.vars = most_important,
     show.overall = F, legend = T, leg.posn = "topleft",
     leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,
     cex.axis = 0.6, line.ylab = 0.9, par.args = list(mgp = c(1.5,
     0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, + 0.3, 0, 0)))
    
    plot(gf, plot.type = "C", imp.vars = most_important,
         show.species = F, common.scale = T, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
       + 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,  0.3, 0, 0)))
    
    plot(gf, plot.type = "P", show.names = T, horizontal = F,
         cex.axis = 1, cex.labels = 0.7, line = 2.5)

    


# RDA on top snps ---------------------------------------------------------
    
    rda_scaled <- flashpcaR::scale2(gen_dat_clim[, top_snps$snp], impute = 2)
    
    # rda_scaled <- flashpcaR::scale2(gen_dat_clim[, snp_col_names], impute = 2)
    
    rda_full <- vegan::rda(rda_scaled ~ latitude + 
                             longitude + elevation +
                             tmax_sum + tmin_winter,
                           data = gen_dat_clim)
    
    rda_full
    summary(rda_full)
    
    anova(rda_full)
    
    anova(rda_full, by="axis", permutations = 100)
    
    anova(rda_full, 
          by = "margin",
          permutations = 100)
    
    round(summary(rda_full)$biplot[, c(1,2)], 2)
    
    plot(rda_full)
    
    
    loadings <- summary(rda_full)$species # Extracts loadings for each SNP (named species here)
    head(loadings)
    
    # Plot histogram of the loadings
    hist(loadings[, 1], xlab = "RDA1", breaks = 30, col = "skyblue2", main = "RDA1")
    
    # Outliers function from Forester et al.
    outliers <- function(loadings, sd_cutoff){
      # find loadings +/-z sd from mean loading
      lims <- mean(loadings) + c(-1, 1) * sd_cutoff * sd(loadings)
      loadings[loadings < lims[1] | loadings > lims[2]] # locus names in these tails
    }
    
    cand_rda1 <- outliers(loadings[, 1], sd_cutoff = 2) 
    cand_rda1
    
    cand_rda2 <- outliers(loadings[, 2], sd_cutoff = 2)
    cand_rda2
    
    outlier_snps_names <- c(names(cand_rda1), names(cand_rda2))
    
    # Extract row indices of outlier snps
    sig_snps <- which(rownames(loadings) %in% outlier_snps_names) 
    
    sig_snps <- which(rownames(loadings) %in% top_snps$snp)
    
    # Color of outlier SNPs
    snp_cols <- rep("grey90",  times = nrow(loadings)) # Set default color
    snp_cols[sig_snps] <- "red"
    
    # Labels of outlier SNPs
    snp_labels <- rep("", times = nrow(loadings))
    snp_labels[sig_snps] <- outlier_snps_names
    
    # Make biplot
    plot(rda_full, type = "n", scaling = 3, xlim = c(-1, 1), ylim = c(-1,1)) # Create empty plot 
    points(rda_full, display = "species", pch = 21,
           col = "black", bg = snp_cols, scaling = 3)
    
    points(rda_full, display = "species", pch = 21, # Overplot sig snps
           col = "black", select = sig_snps, bg = "red", scaling = 3)
    # text(rda_full, display = "species", labels = snp_labels, cex = .75, scaling = 3)
    text(rda_full, scaling = 3, display = "bp", col="black") # Environmental vectors  
    
    
    
    
    # Plot onto map
    
    
    # Get map of California - already loaded in ggplot2 package
    ca_map <- dplyr::filter(ggplot2::map_data("state"), region == "california")
    
    
    # Get loadings for each individual
    ind_loadings <- summary(rda_full)$sites
    head(ind_loadings)
    
    
    # Plot values for first RDA axis
    ggplot(ca_map, aes(x = long, y = lat)) +
      geom_polygon(color = "black", fill = "grey80") +
      geom_point(data = gen_dat_clim, aes(x = longitude, y = latitude, 
                                          fill = ind_loadings[, 1]), # Color by first RDA axis
                 color = "black", pch = 21, size = 5) +
      scale_fill_gradientn(colours = terrain.colors(10), 
                           name = "RDA1") +
      theme_bw() + coord_equal(ratio = 1)
    
    
    
    
    
    

# Simulate data to test out interaction effects in gams -------------------

  dat <- expand.grid(cat = c(1,2,3),
                     x = seq(-3, 3, by = .1))
  dat$cat <- factor(dat$cat)
    
    
  dat$y <- -.5 + -.25*dat$x + -.25*dat$x^2 + rnorm(nrow(dat), sd = .5)
  
  # dat$y[dat$cat == 2] <- -.5 + -.25*dat$x[dat$cat == 2] + .25*dat$x[dat$cat == 2]^2 + rnorm(nrow(dat[dat$cat == 2,]), sd = .5)
  
  # Make cat 1 upside down
  dat$y[dat$cat == 1] <- -.5 + -.25*dat$x[dat$cat == 1] + .25*dat$x[dat$cat == 1]^2 + rnorm(nrow(dat[dat$cat == 1,]), sd = .5)
  dat$y[dat$cat == 2] <- -.5 + -.25*dat$x[dat$cat == 2] + .25*dat$x[dat$cat == 2]^2 + rnorm(nrow(dat[dat$cat == 2,]), sd = .5)
  
  # Make cat 2 a straight line
   dat$y[dat$cat == 2] <- rnorm(nrow(dat[dat$cat == 2,]), sd = .5)
  # dat$y[dat$cat == 3] <- rnorm(nrow(dat[dat$cat == 3,]), sd = .5)
  # 
  ggplot(dat, aes(x, y, col = cat, fill = cat)) + geom_point() + 
    theme_bw() + geom_smooth()
  
  
  # Run model
  gam1 <- gam(y ~ s(x) + cat + s(x, by = cat), data = dat)

  summary(gam1) 
  
  visreg(gam1, xvar = "x", by = "cat")
  
  gam1$coefficients
  
  ## Using s smoother will give a > 0.05 p value when the line is flat
  ## I think m = 1 sets comparison level so that it compares smoothers to each other (https://cran.r-project.org/web/packages/itsadug/vignettes/test.html)
  

