## Notes..

# Re-run models with LGM climate difference and see if significant SNPs overlap

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
  gen_dat_clim <- inner_join(climate_gbs_mom, 
                             dplyr::select(gen_dat, -accession), 
                             by = c("id" = "gbs_name"))
  dim(gen_dat_clim)
  
  climate_gbs_mom$id[!is.na(climate_gbs_mom$accession)]
  
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
    biplot(pca_gen)
    summary(pca_gen)
    screeplot(pca_gen, bstick = TRUE)
    
    pcs <- as.data.frame(pca_gen$x[, 1:10])
    colnames(pcs) <- paste0(colnames(pcs), "_gen")
    
    bglr_gen_scaled <- bind_cols(bglr_gen_scaled, pcs)

  
  # Genetic data as numeric - with mean imputed genotypes
    dat_snp = left_join(dat_gbs_only_scaled,
                        bind_cols(gen_dat_clim[, "accession"],
                                           bglr_gen_scaled),
                        by = "accession")
    dim(dat_snp)
  
  # For using snps as factors
  dat_snp = left_join(dat_gbs_only_scaled,
                    bind_cols(gen_dat_clim[, c("accession", snp_col_names)],
                                       pcs),
                    by = "accession")
  
   
  # Set up factors  
    dat_snp$accession <- factor(dat_snp$accession)
    dat_snp$section_block <- factor(dat_snp$section_block)
    dat_snp$section <- factor(dat_snp$section)
    
    
    # Adding snp_ to column names - useful for gams
    colnames(dat_snp)[which(colnames(dat_snp) %in% snp_col_names)] <- paste0("snp_", colnames(dat_snp[, snp_col_names]))
    
  # Scale PCA axes
    dat_snp$PC1_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC1_gen)) / sd(dat_snp$PC1_gen)
    dat_snp$PC2_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC2_gen)) / sd(dat_snp$PC2_gen)
    dat_snp$PC3_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC3_gen)) / sd(dat_snp$PC3_gen)
    dat_snp$PC4_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC4_gen)) / sd(dat_snp$PC4_gen)
    dat_snp$PC5_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC5_gen)) / sd(dat_snp$PC5_gen)
 
    
  # Scale genotypes  
    scaled_snps_means <- apply(dat_snp[, paste0("snp_", snp_col_names)], MARGIN = 2,
          function(x) mean(x, na.rm = TRUE))
    
    scaled_snps_sds <- apply(dat_snp[, paste0("snp_", snp_col_names)], MARGIN = 2,
                               function(x) sd(x, na.rm = TRUE))
    
    dat_snp[, paste0("snp_", snp_col_names)] <-   apply(dat_snp[, paste0("snp_",
                                             snp_col_names)], MARGIN = 2,
                    function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    
    dat_snp[1:10, 150:160]
    
    # Mean impute
    dat_snp[, paste0("snp_", snp_col_names)] <- apply(dat_snp[, paste0("snp_",
                                                    snp_col_names)], MARGIN = 2,
               function(x) ifelse(is.na(x), yes = 0, no = x))
    
    dat_snp[1:10, 150:160]
    
    summary(dat_snp[, 150:160])
    
    
    
# Set a prediction data frame


 # Prediction frame with a 3 degree increase
  # Set up each climate variable individually
  pred <-  expand.grid(height_2014 = 0,
                  site = "Chico",
                  accession = "1",
                  section_block = "Block1_1",
                  PC1_gen = 0, PC2_gen = 0, PC3_gen = 0)
  
  for(var in climate_vars_dif){
    
    var_unscaled <- paste0(var, "_unscaled")
    
    # Set scaled climate var dif, with a 0 and 3 degree increase as well
    var_df <- expand.grid(accession = "1",
                          var_temp = c(seq(min(dat_all_scaled[ ,var]),
                                               max(dat_all_scaled[ ,var]),
                                               length.out = 50), 
                                           forward_transform(c(0, 3),
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
  # save(dat_snp, snp_col_names, pred,
  #      scaled_snps_means, scaled_snps_sds,
  #      scaled_var_means_gbs_only, scaled_var_sds_gbs_only,
  #   file = paste0("./output/gam_cluster_", Sys.Date(), ".Rdata"))

  
  
  
  
  
# Function for plotting predictions ---------------------------------------
  
  plot_pred <- function(xvar,
                        add_points = TRUE,
                        add_residuals = FALSE,
                        save = FALSE,
                        path = "./",
                        filename = "",
                        main = ""){
    
    # Test to see what type of variable it is - factor or numeric
    type <- ifelse(is.factor(pull(dat_all_scaled, xvar)), "factor", "numeric")
    
    # Make predictions
    v <- visreg(gam_all, xvar = xvar, 
                scale = "response", rug = FALSE, ylab = "5 yr height (cm)", plot = FALSE)
    
    dat_all_scaled_trans <- dat_all_scaled
    
    
    # If xvariable is not a factor, backtranform to original scale for plotting
    if(type == "numeric"){
      # Transform x axis + assign variables for lo and hi CI
      v$fit <- v$fit %>%
        # Back transform X axis
        mutate(!!xvar := back_transform(x = get(xvar),
                                        var = xvar,
                                        means = scaled_var_means_all,
                                        sds = scaled_var_sds_all))
      
      dat_all_scaled_trans <- dat_all_scaled %>%
        # Back transform X axis
        mutate(!!xvar := back_transform(x = get(xvar),
                                        var = xvar,
                                        means = scaled_var_means_all,
                                        sds = scaled_var_sds_all))
    } # End factor if
    
    
    # Add partial residuals
    dat_all_scaled_trans$resid <- v$res$visregRes
    
    ## Main GGPLOT
    gg <- 
      ggplot(dat_all_scaled_trans, aes(x = get(xvar), y = height_2017)) + 
      xlab(xvar) + ylab("5 yr height (cm)") +
      ggtitle(main) +
      theme_bw()
    
    # Add raw points to graph?
    if(add_points == TRUE){
      gg <- gg + 
        geom_jitter(data = dat_all_scaled_trans, aes(x = get(xvar), y = height_2017),
                    size = 0.5, alpha = 0.5)
    }
    
    ## Add partial residuals to graph?
    if(add_residuals == TRUE){
      gg <- gg + 
        geom_point(data = dat_all_scaled_trans, aes(x = get(xvar), y = resid),
                   size = 0.5, alpha = 0.5)
    }
    
    
    # If numeric, add points and lines
    if(type == "numeric"){
      gg <- gg + 
        geom_ribbon(data = v$fit, aes(ymin = visregLwr, ymax = visregUpr), fill = "forestgreen", alpha = 0.75) +
        geom_line(data = v$fit, aes(x = get(xvar), y = visregFit), lwd = 2) 
    }
    
    # If factor, make point_range
    if(type == "factor"){
      gg <- gg + 
        geom_pointrange(data = v$fit, aes(x = get(xvar), y = visregFit, 
                                          ymin = visregLwr, ymax = visregUpr,
                                          fill = get(xvar)),
                        size = 1.5, shape = 21) +
        guides(fill=FALSE) # Remove legend
    }
    
    
    
    # Add verticle line at 0 if a climate transfer function
    if(grepl(pattern = "_dif", x = xvar)){
      gg <- gg + geom_vline(aes(xintercept = 0), lty = 2) 
    }
    
    # General theme and margins
    gg <- gg + theme(plot.margin = (margin(1.5,1.5,1.5,1.5, "cm")),
                     text = element_text(size = 20))
    
    # Printing to file
    if(save){
      ggsave(paste0(path, filename))
    }
    
    print(gg)
    
    return(gg)
    
  } # End plot_pred function
  
  
  # Numerical    
  
  save_flag = FALSE
  path <- "./figs_tables/2018_06_08 Jamie meeting/"
  
  
  
# Run FULL models with all seedlings (no genomic data) --------------------

  # Convert factors
    dat_all_scaled$section <- factor(dat_all_scaled$section)
    dat_all_scaled$section_block <- factor(dat_all_scaled$section_block) 
    dat_all_scaled$accession <- factor(dat_all_scaled$accession)
    dat_all_scaled$site <- factor(dat_all_scaled$site)
    
    var = "tmin_winter_lgm_dif"
    

  for(var in climate_vars_dif){
    
    cat("Working on:", var, "... \n")
    
    form <- formula(paste0("height_2017 ~ section_block + s(height_2014, bs =\"cr\") + s(", var, " , bs = \"cr\") + s(accession, bs = \"re\")"))
    
    
    # With all 5,000+ seedlings
    gam_all <- bam(formula = form,
                   data = dat_all_scaled,
                   discrete = TRUE, 
                   nthreads = 8,
                   method = "fREML", family = "gaussian", 
                   control = list(trace = FALSE))
    
    summary(gam_all)
    
    print(summary(gam_all))
    
    # Plot overall model fit
    test_fit <- dat_all_scaled
    test_fit$pred <- gam_all$fitted.values
    
    
    pdf(paste0("./output/model_visualizations/full models/",var, "_", Sys.Date(), ".pdf" ),
        height = 5, width = 8)

    ggplot(test_fit, aes(x = pred, y = height_2017, bg = section_block)) + geom_point(alpha = 0.75, pch = 21) + theme_bw(15) +
      geom_abline(slope = 1, intercept = 0, lwd = 1.5, col = "forestgreen")

    visreg(gam_all, partial = TRUE)
    visreg(gam_all, partial = FALSE)

    gam.check(gam_all)

   p = plot_pred(xvar = var, add_points = FALSE, add_residuals = FALSE,
              save = FALSE, path = "", filename = "")

  # print(p)

    p =  plot_pred(xvar = var, add_points = FALSE, add_residuals = TRUE,
              save = FALSE, path = "", filename = "")

   # print(p)

    dev.off()

  }  # END LOOP
   
  
# Gam code from cluster for testing -------------------------------------------

  
# Initialize variables
  gam_list <- list()
  x = 1
  
  ## Convert fixed and random parts of formula to actual formula object
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
    snp <- paste0("snp_", snp_col_names[snp_index])
    
    # Choose numbers of knots
    k = length(table(pull(dat_snp, snp)))
    
    # Make sure there's at least 3 levels
    if(length(table(dat_snp[, snp])) < 3){
      cat("Skipping because not at least 3 levels")
      next
    }

    
    # Convert to factor
    # dat_snp[, snp] <- as.factor(pull(dat_snp, snp))
    
  
    # Formula for fixed effects
    
    # Stack overflow on adding m = 1 to by= smooths - https://stats.stackexchange.com/questions/32730/how-to-include-an-interaction-term-in-gam

    
# Run gams  
    cat("Working on: ", snp, "... number: ", x, " ...\n" )    
    
  # For SNPS as continuous  
    fixed_effects_int <- paste0(paste0("height_2017 ~ section_block + s(height_2014, bs=\"cr\") + te(", snp, ", bs=\"cr\", k = ", k, ") +  te(", climate_var_dif,", bs=\"cr\") + ti(", climate_var_dif,", ",snp,", bs=\"cr\", k =", k," )"))
    
    ## Add gen pcs
   # fixed_effects_int <- paste0(fixed_effects_int, "+ s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")")
    
    ## Add MEMs
  #  fixed_effects_int <- paste0(fixed_effects_int,  "+ s(mem1, bs =\"cr\") + s(mem2, bs =\"cr\") + s(mem3, bs =\"cr\") + s(mem4, bs =\"cr\")")
                                       
      # fixed_effects_no_int <- paste0(paste0("height_2017 ~ section_block + s(height_2014, bs=\"cr\") + te(", snp, ", bs=\"cr\", k = ", k, ") +  te(", climate_var_dif,", bs=\"cr\")  + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
    # 
    
    
    
   # Gam with interaction effect
    gam_snp_int = bam(formula = make_formula(fixed_effects_int, random_effects),
                  data = dat_snp, 
                  discrete = TRUE,
                  nthreads = 8,
                  method = "fREML", family = "gaussian", 
                  control = list(trace = FALSE))
    
    # Check for convergence
    if(!gam_snp_int$mgcv.conv){
      cat("Model not converged... trying again \n")
      gam_snp_int = bam(formula = make_formula(fixed_effects_int, random_effects),
                        data = dat_snp, 
                        discrete = TRUE,
                        nthreads = 8,
                        method = "fREML", family = "tw", 
                        control = list(trace = FALSE,
                                       maxit = 500))
      
      if(!gam_snp_int$mgcv.conv){
        cat("Model STILL not converged...\n")
      }
    }

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
    
  # For SNPs as factors  
    # Change accession number if genetic data is missing so we don't get errors in prediction
    # while(any(is.na(dat_snp[as.character(dat_snp$accession) %in% 
    #                         as.character(pred_sub$accession), snp]))){
    #   pred_sub$accession = as.numeric(pred_sub$accession) + 1
    # }
    
    
    ## Visreg plots
    pdf(paste0("./output/model_visualizations/", snp,".pdf"), width = 8, height = 5)

   # With partial residuals
    visreg(gam_snp_int, xvar = climate_var_dif, by = snp,
           overlay = TRUE, partial = TRUE, 
           breaks = c(unique(pred_sub$genotype_scaled)),
           xtrans = function(x) {(x *  scaled_var_sds_gbs_only[climate_var_dif]) + scaled_var_means_gbs_only[climate_var_dif]},
           ylab = "5 year height (cm)")
    
    # Without partial residuals
    visreg(gam_snp_int, xvar = climate_var_dif, by = snp,
           overlay = TRUE, partial = FALSE, rug = FALSE,
           breaks = c(unique(pred_sub$genotype_scaled)),
           xtrans = function(x) {(x *  scaled_var_sds_gbs_only[climate_var_dif]) + scaled_var_means_gbs_only[climate_var_dif]},
           ylab = "5 year height (cm)")
    
    dev.off()

    
    
      
  # Make PREDICTION
  #   gam_preds <-  predict(gam_snp_int, 
  #                         newdata = pred_sub,
  #                         type = "response",
  #                         se = TRUE)
  #   
  #   pred_sub$pred <- gam_preds$fit
  #   pred_sub$se <- gam_preds$se.fit
  #   
  #   pred_sub$genotype <- factor(pred_sub$genotype)
  #   
  #  # Backtransform climate variable
  #    dat_snp_trans <- dat_snp %>%
  #                 dplyr::select(climate_var_dif, snp, height_2017) %>%
  #                    # Back transform X axis
  #                    mutate(!!climate_var_dif := back_transform(x = get(climate_var_dif),
  #                                                    var = climate_var_dif,
  #                                                    means = scaled_var_means_gbs_only,
  #                                                    sds = scaled_var_sds_gbs_only)) %>%
  #                    # Back transform SNP
  #                    mutate(!!snp := back_transform(x = get(snp),
  #                                                               var = snp,
  #                                                               means = scaled_snps_means,
  #                                                               sds = scaled_snps_sds)) %>%
  #                   # Filter down to just 0,1,2 genotypes
  #                   dplyr::filter(get(snp) %in% c(0, 1, 2))
  #    
  #    # Partial residuals
  #    resids <- visreg(gam_snp_int, xvar = climate_var_dif, plot = FALSE)$res
  #    # Back transform X axis
  #    resids <- resids %>%
  #                 mutate(!!climate_var_dif := back_transform(x = get(climate_var_dif),
  #                                               var = climate_var_dif,
  #                                               means = scaled_var_means_gbs_only,
  #                                               sds = scaled_var_sds_gbs_only))
  # 
  # 
  # # Climate transfer functions by genotype   
  #    climate_var_dif_unscaled <- paste0(climate_var_dif, "_unscaled")
  #    
  #   p1 <-  ggplot(pred_sub) + 
  #     geom_ribbon(aes(get(climate_var_dif_unscaled), ymin = pred - 2 * se, ymax = pred + 2 * se, fill = factor(genotype)), 
  #                 alpha = 0.25, show.legend = FALSE)  +
  #     geom_line(data = pred_sub, aes(x = get(climate_var_dif_unscaled), y = pred, col = factor(genotype)), lwd = 1.5) + 
  #     geom_vline(aes(xintercept = 0), lty = 2) +
  #     labs(col = "Genotype") +
  #     ylab("5 yr height (cm)") +
  #     xlab(climate_var_dif) +
  #     geom_point(data = resids, aes(x = get(climate_var_dif), y = visregRes)) +
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
  #   pdf(paste0("./output/model_visualizations/", snp,".pdf"), width = 15, height = 5)  
  #   print(p1 + p2)
  #   dev.off()
  # 
    
  # Make prediction for 3 degree increase in temperature
    for(genotype in c(0,1,2)){
      
      if(!genotype %in% as.numeric(as.character(pred_sub$genotype))){
        gam_snp_int[paste0("height_change_gen_", genotype)] <- NA
        next
      }
      
      pred_sub_temp <- pred_sub[pred_sub$genotype == genotype, ]
      
      height_change <- (pred_sub_temp$pred[pred_sub_temp[, climate_var_dif_unscaled] == 3][1] - 
                          pred_sub_temp$pred[pred_sub_temp[, climate_var_dif_unscaled] == 0][1]) / pred_sub_temp$pred[pred_sub_temp[, climate_var_dif_unscaled] == 3][1] * 100 
      
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
    
    # Converged?
    converged = unlist(lapply(gam_list, function(x) x$mgcv.conv)),
    
    # Deviance explained
    dev_explained = unlist(lapply(gam_list_summary, function(x) x$dev.expl)),
    
    # Rsquared adjusted
    r_sq_adj = unlist(lapply(gam_list_summary, function(x) x$r.sq)),
    
    # Estimated degrees of freedom
    edf = unlist(lapply(gam_list, function(x) sum(x$edf))),
    
    # AIC
    aic = unlist(lapply(gam_list, function(x) x$aic)),
    
    # Delta AIC with model without interaction
    delta_aic = unlist(lapply(gam_list, function(x) x$delta_aic)),
    
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
  
  
  
  

# Read in cluster run output ----------------------------------------------

  path_to_summaries <- "./output/model_summaries_tmax_sum_lgm_dif/"  
  
  sum_df <- ldply(list.files(path_to_summaries, full = TRUE), read_csv)
  
  dim(sum_df)
  sum_df
  
  summary(sum_df)
  
  
  ## Missing snps
  snp_col_names[which(!paste0("snp_", snp_col_names) %in% sum_df$snp)]
  
  
  
  

# Sort and arrange data ---------------------------------------------------

  # How many converged?
  table(sum_df$converged)
  
  # Filter to just those that converged
  sum_df <- sum_df %>%
    dplyr::filter(converged)
  dim(sum_df)
  
  
  
  ## Calculating q values  
      # Vignette - https://bioconductor.org/packages/release/bioc/vignettes/qvalue/inst/doc/qvalue.pdf
      
      qvals <- qvalue(p = sum_df %>%
                        dplyr::pull(p_val_int),
                      fdr.level = 0.05)
      
      summary(qvals)
      
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
    hist(sum_df$p_val_int_adj, breaks = 20)
      
 
  # Calculate average height prediction
    avg_height_pred <- mean(c(pull(sum_df, height_change_gen_0), 
                            pull(sum_df, height_change_gen_1),
                            pull(sum_df, height_change_gen_2)), na.rm = TRUE)
    
    avg_height_pred
    
    avg_height_pred <- 0
  
  # Figuring out which allele is beneficial
    sum_df <- sum_df %>%
        mutate(beneficial_gen_0 = ifelse(height_change_gen_0 > avg_height_pred, 1, 0),
               beneficial_gen_1 = ifelse(height_change_gen_1 > avg_height_pred, 1, 0),
               beneficial_gen_2 = ifelse(height_change_gen_2 > avg_height_pred, 1, 0))
    
    table(sum_df$beneficial_gen_0)
    table(sum_df$beneficial_gen_1)
    table(sum_df$beneficial_gen_2)
    
   # View(sum_df[which(sum_df$beneficial_gen_0 == 1 & sum_df$beneficial_gen_1 == 1 & sum_df$beneficial_gen_2 == 1),])
    
    
  ## Select top snps
    
    # Based on q value
    top_snps = sum_df %>%
    #  dplyr::filter(q_val < 0.05) %>%
    #  dplyr::filter(p_val_int < 0.001) %>%
      dplyr::filter(beneficial_gen_0 == 1 | beneficial_gen_1 == 1 | beneficial_gen_2 == 1) %>%
      dplyr::top_n(-round(nrow(sum_df)*.01), p_val_int) %>% # Take only 1% of snps
     
      dplyr::mutate(snp = gsub("snp_", "", snp)) %>%
      dplyr::arrange(q_val)
    
    top_snps
    
    
  ## Convert top SNPs genotypes to precense absence
    
    gen_dat_ben <- gen_dat_clim
    
    for(snp in top_snps$snp){
      
      sub <- sum_df[sum_df$snp == paste0("snp_", snp),]
      
      ben0 <- sub$beneficial_gen_0 == 1
      ben1 <- sub$beneficial_gen_1 == 1
      ben2 <- sub$beneficial_gen_2 == 1
    
      gen_dat_clim[, snp]
      
      for(row in 1:nrow(gen_dat_clim)){}
      
     gen_dat_ben[, snp] <-  apply(gen_dat_clim[,snp], 1, function(x) {
        if(is.na(x)){
          NA
        } else if(x == 0 & ben0){
          1
        } else if(x == 1 & ben1) {
          1
        } else if(x == 2 & ben2){
          1
        } else {
          0
        }
      })
    }
    
    summary(gen_dat_ben[, top_snps$snp])
    
    
    
    
    
  
  ## SCP code to copy over model plots
    for(snp in paste0("snp_", top_snps$snp)[1:15]){
      system(paste0("scp lukembro@dtn2.hoffman2.idre.ucla.edu:/u/flashscratch/l/lukembro/qlobata_growth/run_291_tmax_sum_dif/model_plots/", snp, "_gg.pdf ./output/model_visualizations/"))
      
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
  
  man_df1 <- data.frame(SNP = paste0("snp_", snp_col_names), 
                       snp_pos)
  
   man_df <- left_join(man_df1, sum_df[, c("snp", "p_val_int")], by = c("SNP" = "snp"))
 
  
  head(man_df)
  
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
  gf_factor <- gen_dat_ben[ ,top_snps$snp]
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

    

# Testing random forest on just beneficial alleles ------------------------

  library(randomForest)    
    
    climate_vars_rf <-  c("tmax_sum", 
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
                           "elevation",
                           "longitude")
    
    pairs.panels(gen_dat_ben[, climate_vars_rf])
    
    
    
    
    
    
    ## Factors
    gf_factor <- gen_dat_ben[ ,top_snps$snp]
    gf_factor <- gf_factor %>%
      replace(is.na(.), 0) %>% # SHould replace with most common case?
      mutate_all(as.factor) 
    
    ## Count number of positives
    pos_thresh <- 50
    gf_factor <- gf_factor[, -which(apply(gf_factor, 2, function(x) sum(x == 1)) < pos_thresh)]
    gf_factor <- gf_factor[, -which(apply(gf_factor, 2, function(x) sum(x == 0)) < pos_thresh)]
    
    dim(gf_factor)
    summary(gf_factor)  
    
    # TO RANDOMIZE GENOTYPES
    # for(col in 1:ncol(gf_factor)){
    #   gf_factor[, col] <- as.factor(sample(0:1, nrow(gf_factor), replace = TRUE))
    # }
    # summary(gf_factor)
    

    
  
    ## Set up rasters
    library(raster)
    
    lobata_range <- readShapePoly("./data/gis/valley_oak_range/qlobata_refined.shp",
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    lobata_range_rough <- readShapePoly("./data/gis/valley_oak_range/querloba.shp",
                                        proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    tmax_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmx1951_1980jja_ave_HST_1513103038/tmx1951_1980jja_ave_HST_1513103038.tif")
    
    dem <- raster("./data/gis/dem/ClimateNA_DEM_cropped.tif")
    
    process_raster <- function(rast){
       rast <- aggregate(rast, fact = 10) # Make smaller
      rast_latlon <- projectRaster(rast, crs = CRS("+proj=longlat +datum=WGS84"))
      plot(rast_latlon)
      return(rast_latlon)
    }
    
    # Aggregate and project tmax_raster
    tmax_rast <- process_raster(tmax_rast)
  
    dem <- resample(dem, tmax_rast) # Match up resolution and extent
    
    compareRaster(tmax_rast, dem) # Make sure rasters are equal

    # Find lat long points
    latlon <- xyFromCell(tmax_rast, 1:ncell(tmax_rast))
    
    rf_rast_df <- data.frame(cell_id = 1:ncell(tmax_rast),
                             longitude = latlon[, 1],
                             latitude = latlon[, 2],
                             tmax_sum = raster::extract(tmax_rast, 1:ncell(tmax_rast)),
                             elevation = raster::extract(dem, 1:ncell(dem)),
                             random = rnorm(1:ncell(dem)))
    rf_rast_df
    
    dim(rf_rast_df)
    
    rf_rast_df_nona <-  rf_rast_df %>%
      filter(!is.na(tmax_sum)) %>%
      filter(!is.na(elevation))
    
    dim(rf_rast_df_nona)
    
  # TO COMPLETELY RANDOMIZE ENVIRONMENTAL VARIABLES
    # gen_dat_ben_randomized <- gen_dat_ben[sample(1:nrow(gen_dat_ben)),]
    # for(var in climate_vars_rf){
    #   gen_dat_ben_randomized[, var] <- rnorm(length(pull(gen_dat_ben_randomized, var)),
    #                                          mean = mean(pull(gen_dat_ben_randomized, var)),
    #                                          sd = sd(pull(gen_dat_ben_randomized, var)))
    # }
    # 
    
  # Begin loop over SNPS  
    stack <- stack()
    oob_error <- NULL
    
    x = 1
    
    for(snp in colnames(gf_factor)){
      
      cat("Working on snp #:", x, "\t", snp, "...\n")
      
      # DO NOT RANDOMIZE CLIMATE VARIABLES
     # X = gen_dat_ben[, climate_vars_rf]
      
      ## RANDOMIZE CLIMATE VARIABLES?
      X = as.data.frame(gen_dat_ben_randomized[, climate_vars_rf])
      
      # Run random forest model
      rf <- randomForest(y = pull(gf_factor, snp),
                         x = X, 
                         importance = TRUE,
                         ntree = 500,
                         mtry = 1)
      
      rf
      
      print(importance(rf))
      print(mean(rf$err.rate[, 1]))

    #  plot(rf)
      
    #  varImpPlot(rf)
      
      # Plot partial plot
      
      # pp_out <- partialPlot(rf, pred.data= X, 
      #                       x.var = "tmax_sum", plot = TRUE)
      # pp_tmax_temp <- data.frame(snp = snp,
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
      
      ## Predict across a raster
      rf_rast_df_nona$pred = predict(rf, rf_rast_df_nona, type = "prob")[, 2] # Second column is positives
      
      rf_rast_df_temp <- left_join(rf_rast_df, rf_rast_df_nona[, c("cell_id", "pred")])

       rast <- tmax_rast
      values(rast) <- rf_rast_df_temp$pred
      
    #  levelplot(rast, margin = FALSE)
      
      stack <- stack(stack, rast)
      
      if(x == 1){
        importance <- data.frame(t(importance(rf)[, 3]))
      }
      
      oob_error[x] <- mean(rf$err.rate[, 1])
      importance <- rbind(importance, importance(rf)[, 3])
      
      x = x + 1
      
    } # End loop over SNPs

    oob_error
    mean(oob_error)
    
    importance
    sort(colMeans(importance), decreasing = TRUE) # Higher is better
     
    stack
    
    mean_stack <-  mean(stack)
    
    levelplot(mean_stack, margin = FALSE)
    
    mean_stack_range <- mask(mean_stack, lobata_range_rough)
    
    levelplot(mean_stack_range, margin = FALSE, contour = FALSE)
    ## ADD ON OBSERVED POINTS FROM SAMPLES OF WHETHER BENEFICIAL ALLELE WAS POSITIVE OR NOT
    ## Make plot 3d based on elevation
    
    
    # Plot partial dependence plots of TMax
    ggplot(pp_tmax, aes(x = tmax_sum, y = tmax_prob,
                        group = snp, alpha = .25)) + geom_line() + theme_bw()
  

  
  
 
 
 

 ## Can maybe overlay a bunch of partial plots of each SNP to show general trends in which direction the gradients are going
    
    
    
    
    
    

# Simulate data to test out interaction effects in gams -------------------

  dat <- expand.grid(cat = c(1,2,3),
                     x = seq(-3, 3, by = .1))
  dat$cat <- factor(dat$cat)
    
    
  dat$y <- -.5 + -.25*dat$x + -.25*dat$x^2 + rnorm(nrow(dat), sd = .5)
  
  # dat$y[dat$cat == 2] <- -.5 + -.25*dat$x[dat$cat == 2] + .25*dat$x[dat$cat == 2]^2 + rnorm(nrow(dat[dat$cat == 2,]), sd = .5)
  
  dat$y[dat$cat == 1] <- -.5 + -.25*dat$x[dat$cat == 1] + .25*dat$x[dat$cat == 1]^2 + rnorm(nrow(dat[dat$cat == 1,]), sd = .5)
  
  ggplot(dat, aes(x, y, col = cat, fill = cat)) + geom_point() + 
    theme_bw() + geom_smooth()
  
  
  gam1 <- gam(y ~ s(x) + cat + s(x, by = cat, m = 1), data = dat)
  summary(gam1) 
  visreg(gam1, xvar = "x", by = "cat")
  

