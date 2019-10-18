

# Load libraries ----------------------------------------------------------

# If running from scratch
# source("./01_clean-process-explore-data.R")
# source("./03_adding-climate-data.R")
#save.image("./output/full growth analysis.Rdata")

# Load R environment that has growth and climate data
load("./output/full growth analysis.Rdata")
dim(dat_all_scaled)

# Load libraries ----------------------------------------------------------
library(tidyverse)
library(googlesheets)
library(maptools)

library(gamm4)
library(visreg)
library(rrBLUP)
library(qvalue) # devtools::install_github("jdstorey/qvalue")
library(MuMIn)
library(vegan)
library(beepr)
library(patchwork) # devtools::install_github("thomasp85/patchwork")
library(rasterVis)
library(pdp)
library(ggrepel)
library(caret)
library(flashpcaR) # devtools::install_github("gabraham/flashpca/flashpcaR")

library(psych)
library(factoextra)

library(foreach)
library(doParallel)
library(gamm4)
library(visreg)
library(MuMIn)
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

# Run FULL models with all seedlings (no genomic data) --------------------

# Convert factors
  dat_all_scaled$section <- factor(dat_all_scaled$section)
  dat_all_scaled$section_block <- factor(dat_all_scaled$section_block, levels = c("IFG_1", "IFG_2", "IFG_3",
                                                                                  "IFG_4", "IFG_5",
                                                                                  "Block1_1", "North annex_2",
                                                                                  "North_3", "North_4", "North_5")) 
  dat_all_scaled$accession <- factor(dat_all_scaled$accession)
  dat_all_scaled$site <- factor(dat_all_scaled$site)
  dat_all_scaled$locality <- factor(dat_all_scaled$locality)
  
  var = "tmax_sum_dif"


# Set formula for gam
  form <- formula(paste0("rgr ~ section_block + s(height_2014, bs =\"cr\") + s(", var, " , bs = \"cr\") + s(PC1_clim, bs =\"cr\") + s(PC2_clim, bs =\"cr\") + s(tmax_sum_dif, by = PC1_clim, bs =\"cr\") + s(tmax_sum_dif, by = PC2_clim, bs =\"cr\") + s(locality, bs = \"re\") + s(accession, bs = \"re\")"))
  
# With all 5,000+ seedlings
  gam_all <- bam(formula = form,
                 data = dat_all_scaled,
                 discrete = FALSE, 
                 nthreads = 8,
                 method = "fREML",
               family = "tw",
               control = list(trace = FALSE))
  
  summary(gam_all)
  
  ## Output data for uploading to figshare
  write_csv(dplyr::select(dat_all_scaled, accession_progeny, accession, locality, section_block,
                          height_2014, tmax_sum_dif, PC1_clim, PC2_clim, rgr), "./output/valley oak growth data - all seedlings.csv")
  
  # visreg(gam_all, partial = TRUE, ylab = "RGR")
  # 
  # visreg(gam_all, partial = FALSE, ylab = "RGR", scale = "response")
  # 
  # gam.check(gam_all)
  
  ## Save gam summary to file
    # sink(file = paste0("./figs_tables/Table S1 - full_gam_model_summary_",
    #                    Sys.Date(), ".txt"))
    # summary(gam_all)
    # anova(gam_all)
    # sink()

  
## Figure S2 - Comparing residual plots to Guassian for Supplementary material
  # gam_all_gaussian <- bam(formula = form,
  #                data = dat_all_scaled,
  #                discrete = TRUE,
  #                nthreads = 8,
  #                method = "fREML",
  #                  family = "gaussian",
  #               # family = "tw",
  #                control = list(trace = FALSE))
  # 
  # summary(gam_all_gaussian)
  # 
  # 
  # par(mfrow = c(1, 2))
  # plot(predict(gam_all, scale = "response"), residuals(gam_all),
  #      xlab = "Linear predictor", ylab = "Residuals", pch = 19,cex = .5, las = 1,
  #      main = "Tweedie error distribution")
  # mtext("(a)", side = 3, line = 1, adj = -0.25, cex = 1.75)
  # plot(predict(gam_all_gaussian), residuals(gam_all_gaussian),
  #      xlab = "Linear predictor", ylab = "Residuals", pch = 19,cex = .5, las = 1,
  #      main = "Gaussian error distribution")
  # mtext("(b)", side = 3, line = 1, adj = -.25, cex = 1.75)
  # 
  # 
  # dev.copy(png, filename = paste0("./figs_tables/Figure S2 - residual plots_",
  #                      Sys.Date(), ".png"),
  #          res = 300, width = 2400, height = 1200)
  # dev.off()

  
  
  # Predicting future increases in temp -------------------------------------
  
  ## Calculate predicted increases in temperature
  future <- read_csv("./data/cleaned_data/maternal tree FUTURE climate data BCM 2018_06_18.csv")
  climate_garden_mom2 <- left_join(climate_garden_mom, future, by = c("accession" = "Accession"))
  
  # High emissions scenario  
  ccsm4_85 <-  climate_garden_mom2$tmax_sum2070_2099jja_ave_CCSM4_rcp85_1529358915 - 
    climate_garden_mom2$tmax_sum
  cnrm_85 <-  climate_garden_mom2$tmax_sum2070_2099jja_ave_CNRM_rcp85_1529358976 - 
    climate_garden_mom2$tmax_sum
  ipsl_85 <-  climate_garden_mom2$tmax_sum2070_2099jja_ave_IPSL_rcp85_1529358959 - 
    climate_garden_mom2$tmax_sum
  miroc_85 <-  climate_garden_mom2$tmax_sum2070_2099jja_ave_MIROC_rcp85_1529357403 - 
    climate_garden_mom2$tmax_sum
  fgoals_85 <- climate_garden_mom2$tmax_sum2070_2099jja_ave_Fgoals_rcp85_1530222451 - 
    climate_garden_mom2$tmax_sum
  
  # RCP 45 scenarios
  miroc_45 <- climate_garden_mom2$tmax_sum2070_2099jja_ave_MIROC_rcp45_1530222488 - 
    climate_garden_mom2$tmax_sum
  mpi_45 <- climate_garden_mom2$tmax_sum2070_2099jja_ave_MPI_rcp45_1530222494 - 
    climate_garden_mom2$tmax_sum
  
  # RCP 26 scenarios
  giss_26 <- climate_garden_mom2$tmax_sum2070_2099jja_ave_GISS_rcp26_1531184801 - 
    climate_garden_mom2$tmax_sum
  miroc5_26 <- climate_garden_mom2$tmax_sum2070_2099jja_ave_MIROC5_rcp26_1531184789 - 
    climate_garden_mom2$tmax_sum
  mri_26 <- climate_garden_mom2$tmax_sum2070_2099jja_ave_MRI_rcp26_1531184812 - 
    climate_garden_mom2$tmax_sum
  
  future_df <-  as.data.frame(cbind(ccsm4_85,
                                    cnrm_85, 
                                    ipsl_85, 
                                    miroc_85,
                                    fgoals_85,
                                    miroc_45,
                                    mpi_45,
                                    giss_26,
                                    miroc5_26,
                                    mri_26))
  
  # pairs.panels(future_df)
  
  ## Density plots of different scenarios for provenances
  ## Miroc much higher than the others
  # future_df %>%
  #   gather(key = "projection", "tmax_sum") %>%
  #   ggplot(., aes(tmax_sum, fill = projection)) + geom_density(alpha = 0.5) + 
  #   theme_bw(15) 
  
  colMeans(future_df)
  
  # Calculate mean increases
  future_85_mean <- mean(apply(dplyr::select(future_df, contains("85")), MARGIN = 1, mean))
  future_45_mean <- mean(apply(dplyr::select(future_df, contains("45")), MARGIN = 1, mean))
  future_26_mean <- mean(apply(dplyr::select(future_df, contains("26")), MARGIN = 1, mean))
  
  
  # Calculate LGM differences in temps
  lgm_mean <- mean(climate_garden_mom$tmax_sum_lgm - climate_garden_mom$tmax_sum)
  
  

  
  

#### Prediction plot of changes in height with tmax transfer  
  
  
  # Make predictions
  v <- visreg(gam_all, xvar = "tmax_sum_dif", 
              scale = "response", rug = FALSE, ylab = "Relative growth rate", plot = FALSE,
              # Set prediction values
              cond = list(section_block = "Block1_1", 
                          locality = "FHL",
                          height_2014 = 0, 
                          PC1_clim = 0,
                          PC2_clim = 0, 
                          accession = 1)) 
  
    # Transform x axis + assign variables for lo and hi CI
    v$fit <- v$fit %>%
      # Back transform X axis
      mutate(tmax_sum_dif = back_transform(x = tmax_sum_dif,
                                      var = "tmax_sum_dif",
                                      means = scaled_var_means_all,
                                      sds = scaled_var_sds_all))
    

## Figure 1 - Plot of decrease in Main GGPLOT ####
    
    
    ## Save version to use as conceptual figure
      # gg <- 
      #   ggplot(v$fit, aes(x = tmax_sum_dif, y = rgr)) +
      #   
      #   geom_vline(aes(xintercept = 0), lty = 2, size = .5) +
      #   scale_x_continuous(breaks = c(0)) +
      #   scale_y_continuous(breaks = NULL) +
      #   ylab("Growth") +
      #   xlab("Tmax climate distance") +
      #   ylim(c(.25, .4)) +
      #   theme_bw(10) + 
      #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      #         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
      #         plot.margin = (margin(1.5,1.5,1.5,1.5, "cm")),
      #         axis.text.y=element_blank(),
      #         axis.ticks.y=element_blank()) +
      #   NULL
      # 
      # gg
      # 
    
    # SAVE TO PDF AND EDIT IN KEYNOTE, etc
    # ggsave(filename = paste0("./figs_tables/fig1/Figure 1 - conceptual diagram ", Sys.Date(), ".pdf"),
    #        gg,
    #        units = "cm", width = 11, height = 8)
    
    
  ## Main plot with results  
   gg <- 
    ggplot(v$fit, aes(x = tmax_sum_dif, y = rgr)) +
    
    geom_vline(aes(xintercept = 0), lty = 2, size = .5) +
     geom_vline(xintercept = lgm_mean, size = .5) + 
    geom_vline(xintercept = future_85_mean, size = .5) +
    geom_ribbon(data = v$fit, aes(ymin = visregLwr, ymax = visregUpr), 
                  fill = "forestgreen", alpha = 0.25) +
    geom_line(data = v$fit, aes(x = tmax_sum_dif, 
                                y = visregFit), lwd = 1,
                                col = "forestgreen") +
    scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5)) +
  #  scale_y_continuous(breaks = seq(.55, .85, .1), limits = c(.55, .85)) +
    ylab(expression(Relative~growth~rate~(cm~cm^-1~yr^-1))) +
  #  ylab("5 yr height (cm)") +
    xlab("Tmax transfer distance (°C) ") +
      # ylim(c(.25, .4)) +
    theme_bw(10) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.margin = (margin(1.5, 1.5, 1.5, 1.5, "cm"))) +
    NULL
  
  gg
  
  
 ## SAVE TO PDF AND EDIT IN KEYNOTE, etc
  # ggsave(filename = paste0("./figs_tables/fig1/Figure 1 - transfer function ", Sys.Date(), ".pdf"),
  #        gg,
  #        units = "cm", width = 11, height = 8)

  # 
  # ggsave(filename = paste0("./figs_tables/Figure S1 - transfer function IFG", Sys.Date(), ".pdf"),
  #        gg,
  #        units = "cm", width = 11, height = 8)
  
  
  
## Plotting site and block effects
  
  
  visreg(gam_all, xvar = "section_block", partial = F, scale = "response", rug = FALSE, 
         ylab = expression(Relative~growth~rate~(cm~cm^-1~yr^-1)), main = "Block effects")
  
  
  
  ## Height
    vis_height <- visreg(gam_all, xvar = "height_2014", partial = F, scale = "response")
  
    vis_height$fit$height_2014 <- back_transform(x =  vis_height$fit$height_2014, 
                                                 var = "height_2014",
                                                 means = scaled_var_means_all,
                                                 sds = scaled_var_sds_all)
  
    gg_height <- ggplot(vis_height$fit, aes(x = height_2014, y = rgr)) +
      
      geom_ribbon(data = vis_height$fit, aes(ymin = visregLwr, ymax = visregUpr), 
                  fill = "black", alpha = 0.25) +
      geom_line(data = vis_height$fit, aes(x = height_2014, 
                                  y = visregFit), lwd = 1,
                col = "black") +
      ylab(expression(Relative~growth~rate~(cm~cm^-1~yr^-1))) +
      xlab("Initial height (cm)") +
      ggtitle("(a) Height effect") + 
      theme_bw(10) + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          #  plot.margin = (margin(1.5, 1.5, 1.5, 1.5, "cm")) 
            ) +
      NULL
    
    gg_height
  
  
  
  ## Blocks
  vis_block <- visreg(gam_all, xvar = "section_block", partial = F, scale = "response")
  
  gg_block <- ggplot(vis_block$fit, aes(x = section_block, y = rgr)) +
    geom_errorbar(data = vis_block$fit, 
                  aes(ymin = visregLwr, ymax = visregUpr), 
                  col = "black", width = .2) +
    geom_point(data = vis_block$fit, aes(x = section_block, 
                                         y = visregFit),
               cex = 2.5,
               col = "black") +
    ylab(expression(Relative~growth~rate~(cm~cm^-1~yr^-1))) +
    geom_vline(xintercept = 5.5, lty = 2) +
    xlab("Block") +
    ggtitle("(b) Block effect") + 
    theme_bw(10) + 
    annotate("text", x = 3, y = .8, label = "Placerville, CA") + 
    annotate("text", x = 8, y = .8, label = "Chico, CA") + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
         # plot.margin = (margin(1.5, 1.5, 1.5, 1.5, "cm")),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
  
  gg_block
  
  
  ## Locality
    vis_locality <- visreg(gam_all, xvar = "locality", partial = F, scale = "response")
    
    gg_locality <- ggplot(vis_locality$fit, aes(x = locality, y = rgr)) +
      geom_errorbar(data = vis_locality$fit, 
                    aes(ymin = visregLwr, ymax = visregUpr), 
                    col = "black") +
      geom_point(data = vis_locality$fit, aes(x = locality, 
                                           y = visregFit),
                 cex = 2.5,
                 col = "black") +
      ylab(expression(Relative~growth~rate~(cm~cm^-1~yr^-1))) +
      xlab("Locality") +
      ggtitle("(c) Locality effect") + 
      theme_bw(10) + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
           # plot.margin = (margin(1.5, 1.5, 1.5, 1.5, "cm")),
            axis.text.x = element_text(angle = 90, hjust = 1, size =  6)) +
      NULL
    
    gg_locality
    
    
    ## Accession
    # vis_accession <- visreg(gam_all, xvar = "accession", partial = F, scale = "response")
    # 
    # gg_accession <- ggplot(vis_accession$fit, aes(x = accession, y = rgr)) +
    #   geom_errorbar(data = vis_accession$fit, 
    #                 aes(ymin = visregLwr, ymax = visregUpr), 
    #                 col = "black") +
    #   geom_point(data = vis_accession$fit, aes(x = accession, 
    #                                           y = visregFit),
    #              cex = 2.5,
    #              col = "black") +
    #   ylab(expression(Relative~growth~rate~(cm~cm^-1~yr^-1))) +
    #   xlab("accession") +
    #   ggtitle("accession effect") + 
    #   theme_bw(10) + 
    #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    #         plot.margin = (margin(1.5, 1.5, 1.5, 1.5, "cm")),
    #         axis.text.x = element_text(angle = 90, hjust = 1)) +
    #   NULL
    # 
    # gg_accession
    
    library(patchwork)
    gg <- gg_height + gg_block + gg_locality + plot_layout(ncol = 1)
    
    ggsave(filename = paste0("./figs_tables/Figure S3 - rgr plots", Sys.Date(), ".pdf"),
           gg,
           units = "cm", width = 17.5, height = 20)
    
    
    

# Estimating heritability of RGR ------------------------------------------

    
    formula(paste0("rgr ~ section_block + s(height_2014, bs =\"cr\") + s(locality, bs = \"re\") + s(accession, bs = \"re\")"))
    
    # With all 5,000+ seedlings
    gam_heritability <- bam(formula =  formula(paste0("rgr ~ section_block + s(height_2014, bs =\"cr\")")),
                   data = dat_all_scaled,
                   discrete = TRUE, 
                   nthreads = 8,
                   method = "fREML", 
                   family = "tw",
                   control = list(trace = FALSE))
    
    hist(residuals(gam_heritability))
    
   lmer_out <-  lmer(residuals(gam_heritability) ~ (1 | accession), 
         data = dat_all_scaled)
   
   var_df <- as.data.frame(VarCorr(lmer_out))
   
   var_va <- var_df$vcov[var_df$grp == "accession"] * (1/(2 * 0.125))
   
   h2 <- var_va / var(residuals(gam_heritability))
   
   h2
  
  
# Predicting changes in height based on degree increase -------------------
  
  # Set degree increases to compare for no change - medium & high emissions & minimum transfer distance
  degrees <- c(0, future_26_mean, future_85_mean, -4.3)
  
  newdata <- data.frame(section_block = v$fit$section_block[1],
                        height_2014 = v$fit$height_2014[1],
                        accession = v$fit$accession[1],
                        locality = v$fit$locality[1],
                        PC1_clim = v$fit$PC1[1],
                        PC2_clim = v$fit$PC2[1],
                        tmax_sum_dif_unscaled = degrees,
                        tmax_sum_dif = forward_transform(x = degrees,
                                                         var = "tmax_sum_dif",
                                                         means = scaled_var_means_all,
                                                         sds = scaled_var_sds_all))
  
  newdata
  
  newdata$pred <- predict(gam_all, newdata = newdata, se.fit = TRUE, type = "response")$fit
  newdata$se <- predict(gam_all, newdata = newdata, se.fit = TRUE, type = "response")$se.fit
  
  newdata$lwr <- newdata$pred - newdata$se * 1.96
  newdata$upr <- newdata$pred + newdata$se * 1.96
  
  newdata
  
  ## % change medium emissions scenario
  
  ## RCP 2.6
  (newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[2]] - 
      newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[1]]) / newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[1]] * 100
  
  
  ## RCP 8.5 change in 5 degree increase
  
  ## Mean
  (newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[3]] - 
      newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[1]]) / newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[1]] * 100
  
  ## Comparing LGM to 
  (newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[1]] - 
      newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[4]]) / newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[4]] * 100
  
  
  
  ## Calculate how much loss in height based on % decrease in RGR
  calc_height <- function(h0, mod, degree_increase, n_years){
   
    
    h1 <- h0
    cat("Height in year zero:", round(h1, 2), " cm with a ", round(degree_increase, 2), 
        " deg_increase: ... \n")
    
    for(year in 1:n_years){
      # Set prediction df
      newdata2 <- data.frame(#section_block = v$fit$section_block[1],
                            section_block = "Block1_1",
                            height_2014 = forward_transform(x = h1,
                                                            var = "height_2014",
                                                            means = scaled_var_means_all,
                                                            sds = scaled_var_sds_all),
                            accession = v$fit$accession[1],
                            locality = v$fit$locality[1],
                            PC1_clim = v$fit$PC1[1],
                            PC2_clim = v$fit$PC2[1],
                            tmax_sum_dif_unscaled = degree_increase,
                            tmax_sum_dif = forward_transform(x = degree_increase,
                                                             var = "tmax_sum_dif",
                                                             means = scaled_var_means_all,
                                                             sds = scaled_var_sds_all))
      
      newdata2
      
      rgr <- predict(gam_all, newdata = newdata2, type = "response")
      print(rgr)
      h1 <- exp(rgr + log(h1))
      
      cat("Height in year:", round(h1, 2), " : ", year, " with deg_increase: ", degree_increase, "... \n")
      
    } # End year loop
    
    return(h1)
  }
  

  h0 = 50

  # Compare these two outputs to determine how much height is lost due to ## % reduction in growth rates under warming climate scenario
  calc_height(h0 = h0, mod = gam_all, 
              degree_increase = 0, 
              n_years = 3)
  
  calc_height(h0 = h0, mod = gam_all, 
              degree_increase = future_85_mean, 
              n_years = 3)
  