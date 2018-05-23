

# Load libraries ----------------------------------------------------------

  # devtools::install_github("thomasp85/patchwork")
    library(patchwork) # For making multi panel ggplots
    library(scales) # For transforming x axis



# Generating own plot with predict

  plot_pred <- function(climate_var_dif, 
                        var2 = NULL, 
                        mod, 
                        base_data,
                        means,
                        sds){
    # Set number of predictions in interval
      length.out = 50
    
    # Set baseline levels
      pred_dat <-  tibble(height_2014 = 0,
                             site = "Chico",
                             accession = "1",
                             section_block = "Block1_1")
    
    ## New data frame with sequence of variable in question
      sequence <- data.frame(accession = pred_dat$accession[1],
                 climate_var_dif = seq(from = min(base_data[, climate_var_dif]),
                                       to = max(base_data[, climate_var_dif]),
                                       length.out = length.out))
      
      # Rename climate_var_dif to actual name of variable
      sequence <- sequence %>% mutate(!!climate_var_dif := climate_var_dif,
                                      climate_var_dif = NULL)
      
      # Join to prediction data frame
      pred_dat <- left_join(pred_dat, sequence)
      pred_dat
      
      
    # If 2d plot - do again for second variable
      if(!is.null(var2)){
        sequence2 <- data.frame(accession = pred_dat$accession[1],
                               var2 = seq(from = min(base_data[, var2]),
                                                     to = max(base_data[, var2]),
                                                     length.out = length.out))
        # Rename var2 to actual name
        sequence2 <- sequence2 %>% mutate(!!var2 := var2,
                                        var2 = NULL)
        
        # Join to main dataframe
        pred_dat <- left_join(pred_dat, sequence2)
        pred_dat
        
      }
     
    
    
    # Make predictions with SE
      predictions <- predict(mod, newdata = pred_dat, 
                             type = "response",
                             se.fit = TRUE)
    
      pred_dat$pred <- predictions$fit
      pred_dat$se <- predictions$se.fit
    
    # Function to back transform variables
    back_transform <- function(x, var, means, sds){
      (x * sds[var]) + means[var]
    }
    
    # Transform x axis + assign variables for lo and hi CI
    pred_dat <- pred_dat %>%
      # Back transform X axis
      mutate(!!climate_var_dif := back_transform(x = get(climate_var_dif),
                                                 var = climate_var_dif,
                                                 means = means,
                                                 sds = sds)) %>%
                            mutate(pred_hi = pred + 2 * se, 
                                   pred_lo = pred - 2 * se)
    
    # Create ggplot of predictions
     #  p = pred_dat %>%
     #    # Back transform X axis
     #    mutate(!!climate_var_dif := back_transform(x = get(climate_var_dif), 
     #                                            var = climate_var_dif,
     #                                            means = scaled_var_means_all,
     #                                            sds = scaled_var_sds_all)) %>%
     #  ggplot(., aes(x = get(climate_var_dif), y = pred)) +
     # geom_ribbon(aes(ymin = pred - 2 * se, ymax = pred + 2 * se),
     #              fill = rgb(35, 139, 35, alpha = 50, maxColorValue = 255)) +
     #   geom_line(lwd = 2, col = rgb(35, 139, 35, maxColorValue = 255)) +
     #  geom_vline(xintercept = 0, lty = 2) +
     #   ylab("5 yr height (cm)") + xlab(climate_var_dif) +
     # theme_bw()
     #  
     #  print(p)
     #  
     #  return(p)
    
    
    ## 2d basic plotting
    
    if(is.null(var2)){
      
      par(mar=c(5.1,4.1,2.1,2.1))
        
        # Base empty plot
        plot(x = pull(pred_dat, climate_var_dif), y = pred_dat$pred,
             xlab = climate_var_dif, ylab = "5 yr height (cm)",
             type = "n", las = 1, ylim = c(min(pred_dat$pred_lo), max(pred_dat$pred_hi)))
        
        # 0 line for no transfer
        abline(v = 0, lty = 2)
        
        # 95% CI
        polygon(c(pull(pred_dat, climate_var_dif), rev(pull(pred_dat, climate_var_dif))),
                c(pull(pred_dat, pred_lo), rev(pull(pred_dat, pred_hi))), 
                col=rgb(35, 139, 35, alpha = 50, maxColorValue = 255), border=NA)
        
        # Main prediction line
        lines(x = pull(pred_dat, climate_var_dif), y = pred_dat$pred, 
             col = rgb(35, 139, 35, maxColorValue = 255), 
             lwd = 2)
        
    }
      
      
   
   ## 3d perspective plot
    if(!is.null(var2)){
    
      par(mar=c(1,1,1,1)+0.1)
      
    # Make list of unique increasing x and y values for axes   
      pred_dat2 <- pred_dat %>%
        dplyr::select(climate_var_dif) %>%
        unique() %>%
        bind_cols(., pred_dat %>%
                    dplyr::select(var2) %>%
                    unique())
  
      # Munge into persp format - an x*y matrix with z values
        persp_format <- pred_dat %>%
                      dplyr::select(climate_var_dif, var2, pred) %>%
                      spread(key = var2, value = pred)
        
        persp_format
        
        # Remove first column - climate transfer
        persp_format <- as.matrix(persp_format[, -1])
      
      
      # Create a function interpolating colors in the range of specified colors
        # From persp help function
        
        alpha_surface = 255
        
        color_ramp <- colorRampPalette(c(rgb(43,130,189, alpha_surface, maxColorValue = 255),
                                         rgb(247,247,247, alpha_surface, maxColorValue = 255),
                                         rgb(35, 139, 35, alpha_surface, maxColorValue = 255)),
                                       alpha = TRUE)
       
       # Generate the desired number of colors from this palette
        nbcol <- 100
        color <- color_ramp(nbcol)
        # Compute the z-value at the facet centres
        zfacet <- persp_format[-1, -1] + 
                  persp_format[-1, -ncol(persp_format)] + 
                  persp_format[-nrow(persp_format), -1] +
                  persp_format[-nrow(persp_format), -ncol(persp_format)]
        # Recode facet z-values into color indices
        facetcol <- cut(zfacet, nbcol)
      
     persp(x = pull(pred_dat2, climate_var_dif),
            y = pull(pred_dat2, var2),
            z = persp_format,
            zlab = "\n\n5 year height (cm)", 
            xlab = paste0("\n", climate_var_dif),
            ylab = paste0("\n", var2),
            theta = 30, phi = 15,
            d = 2,  
            ticktype = "detailed",
            border = rgb(0,0,0, alpha = 50, maxColorValue = 255),
          #  col = rgb(35, 139, 35, alpha = 175,  maxColorValue = 255),
          col = color[facetcol]) -> p
     
     base_data2 <- base_data %>%
       mutate(!!climate_var_dif := back_transform(x = get(climate_var_dif),
                                                var = climate_var_dif,
                                                means = means,
                                                sds = sds))

     points(trans3d(x = pull(base_data2, climate_var_dif),
                     y = pull(base_data2, var2),
                     z = base_data$height_2017, pmat = p), pch = 19, cex = .25)
      
    } # End 3d plot IF
  
  } # End function



# Make plots --------------------------------------------------------------


### Figure 2 - 
  pdf(file = paste0("./figs_tables/Figure2_", Sys.Date(), ".pdf"))
  
  par.default <- par(no.readonly=TRUE)
  par(mfrow = c(3,2))

# Tmax_summer_dif    
  plot_pred(climate_var_dif = "tmax_sum_dif", 
            mod = gam_mods$tmax_sum$gam_clim_dif_all,
            base_data = dat_all_scaled,
            means = scaled_var_means_all, sds = scaled_var_sds_all)
  
  # 3d
  plot_pred(climate_var_dif = "tmax_sum_dif", 
            var2 = "bv.tmax_sum.kin",
            mod = gam_mods$tmax_sum$gam_kin_int,
            base_data = dat_bv,
            means = scaled_var_means_gbs_only, sds = scaled_var_sds_gbs_only)
  
  plot_pred(climate_var_dif = "tmax_sum_dif", 
            var2 = "bv.tmax_sum.kin.marker100",
            mod = gam_mods$tmax_sum$gam_kin_marker100_int,
            base_data = dat_bv,
            means = scaled_var_means_gbs_only, sds = scaled_var_sds_gbs_only)
  
# Tmin winter dif  
  plot_pred(climate_var_dif = "tmin_winter_dif", 
            mod = gam_mods$tmin_winter$gam_clim_dif_all,
            base_data = dat_all_scaled,
            means = scaled_var_means_all, sds = scaled_var_sds_all)
  # 3d
  plot_pred(climate_var_dif = "tmin_winter_dif", 
            var2 = "bv.tmin_winter.kin",
            mod = gam_mods$tmin_winter$gam_kin_int,
            base_data = dat_bv,
            means = scaled_var_means_gbs_only, sds = scaled_var_sds_gbs_only)
  
  
# Tmin winter dif  
  plot_pred(climate_var_dif = "bioclim_04_dif", 
            mod = gam_mods$bioclim_04$gam_clim_dif_all,
            base_data = dat_all_scaled,
            means = scaled_var_means_all, sds = scaled_var_sds_all)
  # 3d
  plot_pred(climate_var_dif = "bioclim_04_dif", 
            var2 = "bv.bioclim_04.kin",
            mod = gam_mods$bioclim_04$gam_kin_int,
            base_data = dat_bv,
            means = scaled_var_means_gbs_only, sds = scaled_var_sds_gbs_only)  

  
  par(par.default)

  dev.off()
  
  
 test =  bam(height_2017 ~ site + section_block + s(height_2014) + te(tmax_sum_dif) + 
        te(bv.PC1.kin) + ti(tmax_sum_dif, bv.PC1.kin) + 
        s(accession, bs = "re"), data = dat_bv, 
        nthreads = 8,
        method = "fREML", family = "tw")
 
 summary(test)
 
 visreg2d(test, xvar = "tmax_sum_dif", yvar = "bv.PC1.kin", scale = "response", plot.type = "persp")
 
 AICc(gam_mods$tmax_sum$gam_clim_dif)
 AICc(test)
 
 

# Testing out individual SNP models ---------------------------------------------------------

 
 dat_bv2 <- dat_bv
 dat_bv2$accession <- as.numeric(as.character(dat_bv$accession))
 
 dat_snp = left_join(dat_bv2, bind_cols(gen_dat_clim[, "accession"],
                                        bglr_gen_scaled), 
                     by = "accession")
 
 dat_snp$accession <- factor(dat_snp$accession)
 

 counts <- apply(dat_snp[, snp_col_names], 2, function(x) length(table(x)))
 
 snp_col_names_sub <- snp_col_names[counts >= 3 ]
 
 save(dat_snp, snp_col_names, file = "./output/gam_cluster.Rdata")
 
 length(snp_col_names)
 
 
 
 
 # Initialize variables
 gam_list <- list()
 x = 1
 
 
 ## Set up formulas for gam models
 climate_var_dif = "tmax_sum_dif"
 
 
 ## Convert fixed and random parts of formula to actual formula objetc
 # No idea why we need to call formula twice, but it works
 make_formula <- function(fixed, random){
   return(formula(formula(enquote(paste0(fixed, random)))))
 }
 
 # Formula for fixed effects
 fixed_effects <- paste0(paste0("height_2017 ~ site + s(height_2014) + te(", 
                                climate_var_dif,") + te(snp_dbl, k = 3) + ti(",
                                climate_var_dif, ", snp_dbl, k = 3)"))
 
 # Formula for random effects
 random_effects <-  '+ s(accession, bs = "re") + s(section_block, bs = "re")'
 
 task_id = 10; 
 interval = 10
 
 # Loop through snps by index based on task id
 for(snp_index in task_id:(task_id + interval - 1)){
   
   # For timing loops
   start_time <- Sys.time()
   
   # Choose snp
   snp <- snp_col_names_sub[snp_index]
   
   cat("Working on: ", snp, "...\n" )
   
   #dat_snp$snp_factor <- factor(pull(dat_snp, snp))
   dat_snp$snp_dbl <- pull(dat_snp, snp)
   
   gam_snp = bam(formula = make_formula(fixed_effects, random_effects),
                 data = dat_snp, 
                 discrete = TRUE,
                 nthreads = 8,
                 method = "fREML", family = "tw")
   
   # Save into list
   gam_list[[x]] <- gam_snp
   names(gam_list)[x] <- snp
   
   x = x + 1
   
   end_time <- Sys.time()
   
   cat("This loop took: ")
   print(end_time - start_time)
   
   # Save individual .Rdata file for each model
   cat("Saving gam model to file... \n")
   
   gam_name <- paste0("gam_", climate_var_dif, "_", "snp")
   assign(x = gam_name, value = gam_snp) # Assign model the name
   
   save(list = gam_name, file = paste0("./", gam_name, ".Rdata"))
   
 }
 
 # Model summaries
 
 # SNP name
 
 
 # Run summary for each model
 gam_list_summary <- lapply(gam_list, summary)
 
 gam_mod_summary_df <- data.frame(
   
            # climate var
            climate_var = climate_var_dif,
   
            # Snp name
            snp = names(gam_list),
            # Deviance explained
            dev_explained = unlist(lapply(gam_list_summary, function(x) x$dev.expl)),
            
            # Estimated degrees of freedom
            edf = unlist(lapply(gam_list, function(x) sum(x$edf))),
            
            # AIC
            aic = unlist(lapply(gam_list, function(x) x$aic)),
            
            # P value of interaction term
            p_val_int = unlist(lapply(gam_list_summary, function(x) x$s.table[paste0("ti(", climate_var_dif, ",snp_dbl)"), "p-value"]))
            
            )
 
 gam_mod_summary_df
 
 write.csv(gam_mod_summary_df, file = "./")

 
 
 
 
 
 snp_top <- "1_243802"
 
 summary(sig_gams[[snp_top]])
 
 plot(dat_snp$tmax_sum_dif, pull(dat_snp, snp_top))
 
 visreg(sig_gams[[snp_top]], xvar = "tmax_sum_dif", by = "snp_dbl", scale = "response")
 

 visreg2d(sig_gams[[snp_top]], xvar = "tmax_sum_dif", yvar = "snp_dbl", scale = "response", plot.type = "persp", theta = 45, phi = 15)
 
 plot(dat_snp$tmax_sum_dif, pull(dat_snp, snp_top))
 
 pairs.panels(dplyr::select(dat_snp, snp_top, climate_vars))
 
 AIC(sig_gams[[snp_top]])
 AIC(gam_base)
 

 gam_base = bam(height_2017 ~ 
                 site +
                 + s(height_2014, k = 5) 
               + s(tmax_sum_dif) 
               + s(section_block, bs = "re")
               + s(accession, bs = "re"), 
               data = dat_snp, 
               discrete = TRUE, nthreads = 8,
               method = "fREML", family = "tw")
 
 visreg(gam_base, scale = "response")
 AIC(gam_base)
 
 
