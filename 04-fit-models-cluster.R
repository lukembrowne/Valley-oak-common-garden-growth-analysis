## Notes..

# Switch to gaussian to better display partial residual plots - can then show that log transforming doesn't change anything
# Re-run models with LGM climate difference and see if significant SNPs overlap

# Load libraries ----------------------------------------------------------

# install_github("jdstorey/qvalue")
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

  
  
  
  
# Run FULL models with all seedlings (no genomic data) --------------------

  # Convert factors
    dat_all_scaled$section <- factor(dat_all_scaled$section)
    dat_all_scaled$section_block <- factor(dat_all_scaled$section_block) 
    dat_all_scaled$accession <- factor(dat_all_scaled$accession)
    dat_all_scaled$site <- factor(dat_all_scaled$site)
  
  # With all 5,000+ seedlings
  gam_all <- bam(height_2017 ~ section_block + 
                          + s(height_2014, bs= "cr")  
                          + s(tmax_sum_dif, bs="cr")
                          + s(accession, bs = "re"),
                          data = dat_all_scaled,
                          discrete = TRUE, 
                          nthreads = 8,
                          method = "fREML", family = "gaussian", 
                          control = list(trace = TRUE))
  
  summary(gam_all)

# Plot overall model fit
  test_fit <- dat_all_scaled
  test_fit$pred <- gam_all$fitted.values
  
  ggplot(test_fit, aes(x = pred, y = height_2017, bg = section_block)) + geom_point(alpha = 0.75, pch = 21) + theme_bw(15) + 
    geom_abline(slope = 1, intercept = 0, lwd = 1.5, col = "forestgreen") 
  
  visreg(gam_all, partial = TRUE)
  visreg(gam_all, partial = FALSE)
  
  gam.check(gam_all)

  
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
      
  } # End plot_pred function

    
# Numerical    
  
  save_flag = TRUE
  path <- "./figs_tables/2018_06_08 Jamie meeting/"
  
  plot_pred(xvar = "height_2014", add_points = FALSE, save = save_flag, path = path, filename = "height_2014_no_points.pdf")
  plot_pred(xvar = "height_2014", add_points = TRUE, save = save_flag, path = path, filename = "height_2014_points.pdf")
  
  plot_pred(xvar = "tmax_sum_dif", add_points = FALSE, add_residuals = TRUE, main = "Tmax _summer 1950-1980", 
            save = save_flag, path = path, filename = "tmax_sum_dif_no_points.pdf")
  plot_pred(xvar = "tmax_sum_dif", add_points = FALSE, main = "Tmax _summer 1950-1980", 
            save = save_flag, path = path, filename = "tmax_sum_dif_points.pdf")
  
# Historic  
  plot_pred(xvar = "tmax_sum_last1000_dif", add_points = FALSE, add_residuals = TRUE,
            main = "Tmax_summer Last Millennium (1,000 yrs ago)",
            save = save_flag, path = path, filename = "tmax_sum_last1000_dif_no_points.pdf")
  plot_pred(xvar = "tmax_sum_last1000_dif", add_points = FALSE, main = "Tmax_summer Last Millennium (1,000 yrs ago)",
            save = save_flag, path = path, filename = "tmax_sum_last1000_dif_points.pdf")
  
  plot_pred(xvar = "tmax_sum_holo_dif", add_points = FALSE, add_residuals = TRUE, main = "Tmax_summer Mid-Holocene (6,000 yrs ago)",
            save = save_flag, path = path, filename = "tmax_sum_holo_dif_no_points.pdf")
  plot_pred(xvar = "tmax_sum_holo_dif", add_points = TRUE, main = "Tmax_summer Mid-Holocene (6,000 yrs ago)",
            save = save_flag, path = path, filename = "tmax_sum_holo_dif_points.pdf")
  
  plot_pred(xvar = "tmax_sum_lgm_dif", add_points = FALSE, add_residuals = TRUE, 
            main = "Tmax_summer Last Glacial Maximum (21,000 yrs ago)",
            save = save_flag, path = path, filename = "tmax_sum_lgm_dif_no_points.pdf")
  
  plot_pred(xvar = "tmax_sum_lgm_dif", add_points = FALSE, add_residuals = FALSE, 
            main = "Tmax_summer Last Glacial Maximum (21,000 yrs ago)",
            save = save_flag, path = path, filename = "tmax_sum_lgm_dif_points.pdf")
  
  
 # Tmin winter 
  plot_pred(xvar = "tmin_winter_dif", add_points = FALSE, save = save_flag, path = path, filename = "tmin_winter_dif_no_points.pdf")
  plot_pred(xvar = "tmin_winter_dif", add_points = TRUE, save = save_flag, path = path, filename = "tmin_winter_dif_points.pdf")
  
  
  plot_pred(xvar = "tmin_winter_last1000_dif",  add_points = FALSE, add_residuals = TRUE, 
            main = "Tmin_winter Last Millenium (1,000 yrs ago)",
            save = save_flag, path = path, filename = "tmin_winter_last1000_dif_no_points.pdf")
  plot_pred(xvar = "tmin_winter_last1000_dif", add_points = FALSE,  add_residuals = FALSE, 
            main = "Tmin_winter Last Millenium (1,000 yrs ago)",
            save = save_flag, path = path, filename = "tmin_winter_last1000_dif_points.pdf")
  
  plot_pred(xvar = "tmin_winter_holo_dif", add_points = FALSE, 
            save = save_flag, path = path, filename = "tmin_winter_holo_dif_no_points.pdf")
  plot_pred(xvar = "tmin_winter_holo_dif", add_points = TRUE, 
            save = save_flag, path = path, filename = "tmin_winter_holo_dif_points.pdf")
  
  plot_pred(xvar = "tmin_winter_lgm_dif", add_points = FALSE, add_residuals = TRUE, 
            main = "Tmin_winter Last Glacial Maximum (21,000 yrs ago)",
            save = save_flag, path = path, filename = "tmin_winter_lgm_dif_no_points.pdf")
  plot_pred(xvar = "tmin_winter_lgm_dif", add_points = FALSE,  add_residuals = FALSE, 
            main = "Tmin_winter Last Glacial Maximum (21,000 yrs ago)",
            save = save_flag, path = path, filename = "tmin_winter_lgm_dif_points.pdf")
  
# Factors
  plot_pred(xvar = "section_block", save = save_flag, path = path, filename = "section_block.pdf")  
  plot_pred(xvar = "site", save = save_flag, path = path, filename = "site.pdf")  
  plot_pred(xvar = "accession", add_points = FALSE, save = save_flag, path = path, filename = "accession.pdf")
      
      
  
# Gam code from cluster for testing -------------------------------------------

  
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
    
  # For SNPS as factors    
    
  # For factor interaction  
    # fixed_effects_int <- paste0(paste0("height_2017 ~ site + ", snp, " + s(height_2014, bs=\"cr\") + s(", climate_var_dif,", bs=\"cr\") + s(", climate_var_dif,", by =",snp,", bs=\"cr\", m = 1) + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
    
    # fixed_effects_no_int <- paste0(paste0("height_2017 ~ site + ", snp, " + s(height_2014, bs=\"cr\") + s(", climate_var_dif,", bs=\"cr\")  + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
    
    
    
# Run gams  
    cat("Working on: ", snp, "... number: ", x, " ...\n" )    
    
  # For SNPS as continuous  
    fixed_effects_int <- paste0(paste0("height_2017 ~ section_block + s(height_2014, bs=\"cr\") + te(", snp, ", bs=\"cr\", k = ", k, ") +  te(", climate_var_dif,", bs=\"cr\") + ti(", climate_var_dif,", ",snp,", bs=\"cr\", k =", k," ) + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
    
    fixed_effects_no_int <- paste0(paste0("height_2017 ~ section_block + s(height_2014, bs=\"cr\") + te(", snp, ", bs=\"cr\", k = ", k, ") +  te(", climate_var_dif,", bs=\"cr\")  + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
    
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

  path_to_summaries <- "./output/model_summaries_tmax_sum/"  
  
  sum_df <- ldply(list.files(path_to_summaries, full = TRUE), read_csv)
  
  dim(sum_df)
  sum_df
  
  summary(sum_df)
  
  
  
  

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
      
        qvals$qvalues
        qvals$lfdr
        qvals$pi0 # overall proportion of true null hypotheses
      
      # Sort by lowest q value  
      sum_df %>%
        dplyr::arrange(q_val)
      
  
 
      
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
      dplyr::filter(q_val < 0.05) %>%
    #  dplyr::filter(beneficial_gen_0 == 1 | beneficial_gen_1 == 1 | beneficial_gen_2 == 1) %>%
      #  dplyr::select(snp) %>%
      dplyr::mutate(snp = gsub("snp_", "", snp)) %>%
      dplyr::arrange(q_val)
    
    top_snps
  
  ## SCP code to copy over model plots
    for(snp in paste0("snp_", top_snps$snp)){
      system(paste0("scp lukembro@dtn2.hoffman2.idre.ucla.edu:/u/flashscratch/l/lukembro/qlobata_growth/run_2876655_tmax_sum_dif/model_plots/", snp, ".pdf ./output/model_visualizations/"))
      
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
         chr.den.col=c("darkgreen", "yellow", "red"), file.output = TRUE)
  
## QQplot  
  dev.off()
  CMplot(man_df, plot.type="q", box=TRUE,file.output = TRUE)
  
# With Q values  
   man_df <- left_join(man_df1, sum_df[, c("snp", "q_val")], by = c("SNP" = "snp"))
   dev.off()
   CMplot(man_df, plot.type=c("m"), LOG10=TRUE, threshold= 0.05,
          bin.size=1e6, ylim = c(0, 2),
          chr.den.col=c("darkgreen", "yellow", "red"), file.output = TRUE)
  

  
  

# RDA on top snps ---------------------------------------------------------

  rda_scaled <- flashpcaR::scale2(gen_dat_clim[, top_snps$snp], impute = 2)
  

  
  rda_full <- vegan::rda(rda_scaled ~ latitude + 
                           longitude + elevation + tmax_sum + tmin_winter,
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
  
  # Color of outlier SNPs
  snp_cols <- rep("grey90",  times = nrow(loadings)) # Set default color
  snp_cols[sig_snps] <- "red"
  
  # Labels of outlier SNPs
  snp_labels <- rep("", times = nrow(loadings))
  snp_labels[sig_snps] <- outlier_snps_names
  
  # Make biplot
  plot(rda_full, type = "n", scaling = 3) # Create empty plot 
  points(rda_full, display = "species", pch = 21,
         col = "black", bg = snp_cols, scaling = 3)
  text(rda_full, display = "species", labels = snp_labels, cex = .75, scaling = 3)
  text(rda_full, scaling = 3, display = "bp", col="black") # Environmental vectors  
  
  
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
  
  gf_scaled <- flashpcaR::scale2(gen_dat_clim[, top_snps$snp], impute = 2)
   
    
    library(gradientForest)
    
    climate_vars_gf <- c("tmax_sum", "tmin_winter", "bioclim_04", 
                         "random", "cwd", "ppt", "longitude", "elevation",
                         "latitude")
    
    pairs.panels(gen_dat_clim[, climate_vars_gf])
    
    length(top_snps$snp)
    top_snps$snp
    
    gf <- gradientForest(cbind(gen_dat_clim[, climate_vars_gf],
                                     bglr_gen_scaled ),
                         predictor.vars = climate_vars_gf, 
                       # response.vars = top_snps$snp,
                       response.vars = gsub("snp_", "", sum_df$snp),
                        trace = TRUE)

    gf$result # Rsquared of positive loci
    gf$species.pos.rsq
    
    mean(gf$result)
    
    plot(gf, plot.type = "O")
    
    most_important <- names(importance(gf))
    
    pred <- expand.grid(latitude = seq(35, 41, by = .1))
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
    
    plot(gen_dat_clim$latitude, pull(gen_dat_clim, snp))
    
    
    
    
    
    
    
    
    

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
  

