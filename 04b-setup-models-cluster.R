

# Load libraries ----------------------------------------------------------
  
  source("./01_clean-process-explore-data.R")
  source("./03_adding-climate-data.R")
  
  dim(dat_all_scaled)
  dim(dat_gbs_all_scaled)

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
  
  
  # function to transform from raw to scaled variables
    forward_transform <- function(x, var, means, sds){
      ( x - means[var]) / sds[var]
    }
    
    # Function to back transform variables
    back_transform <- function(x, var, means, sds){
      (x * sds[var]) + means[var]
    }

# Format and save to .Rdata for cluster -----------------------------------

  # Do PCA on kinship matrix  
    
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
  # missing <- gen_dat_clim %>%
  #   dplyr::select(id, accession) %>%
  #   mutate(missing = apply(gen_dat_clim[, snp_col_names], 1, function(x) sum(is.na(x))/length(x)))
  # 
  # missing %>%
  #   arrange(desc(missing))

  # Scale data and impute missing before converting 012 to -101
  bglr_gen_scaled <- flashpcaR::scale2(gen_dat_clim[, snp_col_names], impute = 2)
  bglr_gen_scaled[1:10, 1:10]
  summary(bglr_gen_scaled[, 1:10])
  
  
  ## Calculate kinship matrix
  
  # Need to first convert to -1, 0, 1 matrix
    kin_mat <- as.matrix(gen_dat_clim[, snp_col_names])
    kin_mat[kin_mat == 0] <- -1
    kin_mat[kin_mat == 1] <- 0
    kin_mat[kin_mat == 2] <- 1
    
    kin_mat <- A.mat(kin_mat) # Calculate kinship matrix
  
  
  ## Calculate PCA on kinship matrix
    pca_gen = prcomp(kin_mat, center = TRUE, scale = TRUE) ## PCA on kinship matrix
  
   # biplot(pca_gen)
    summary(pca_gen)
  #  screeplot(pca_gen, bstick = TRUE)
    
  ### Figure S3 - Plotting results of PCA for supplementary material

      # par(mfrow = c(2,1))
      # 
      # # Scale values for plotting in biplot - from: https://stats.stackexchange.com/questions/66926/what-are-the-four-axes-on-pca-biplot
      # scores<-pca_gen$x
      # lam <- pca_gen$sdev[1:3]
      # n <- NROW(scores)
      # lam <- lam * sqrt(n)
      # lam <- lam^1
      # yy<-t(t(pca_gen$rotation[, 1:3]) * lam)
      # xx<-t(t(scores[, 1:3])/lam)
      # 
      # # Pc1 vs pc2
      #   plot(xx[, 1], xx[, 2],
      #        las = 1,
      #        xlab = paste0("PC1: ", summary(pca_gen)$importance[2, 1] * 100, "% variance explained"),
      #        ylab = paste0("PC2: ", summary(pca_gen)$importance[2, 2] * 100, "% variance explained"),
      #        pch = 19,
      #        type = "n")
      #   text(xx[, 1], xx[, 2],
      #        labels = gen_dat_clim$id, cex = .5)
      #   mtext(text = "(a)", side = 3, adj = -.1, padj = -1, cex = 1.25)
      # 
      # 
      #   # Pc1 vs pc3
      #   plot(xx[, 1], xx[, 3],
      #        las = 1,
      #        xlab = paste0("PC1: ", summary(pca_gen)$importance[2, 1] * 100, "% variance explained"),
      #        ylab = paste0("PC3: ", summary(pca_gen)$importance[2, 3] * 100, "% variance explained"),
      #        pch = 19,
      #        type = "n")
      #   text(xx[, 1], xx[, 3],
      #        labels = gen_dat_clim$id, cex = .5)
      #   mtext(text = "(b)", side = 3, adj = -.1, padj = -1, cex = 1.25)
      # 
      #   dev.copy(png, filename = paste0("./figs_tables/Figure S3 - PCA plot_",
      #                                   Sys.Date(), ".png"),
      #            res = 300, width = 1800, height =2400)
      #   dev.off()

    pcs <- as.data.frame(pca_gen$x[, 1:10])
    colnames(pcs) <- paste0(colnames(pcs), "_gen")
    
    bglr_gen_scaled <- bind_cols(bglr_gen_scaled, pcs)

    
 
    
 # Set up dataframe with SNP data, height data, and climate data
    dat_snp_all = left_join(dat_gbs_all_scaled,
                      bind_cols(gen_dat_clim[, c("accession", snp_col_names)],
                                         pcs),  by = "accession")
    
  # Set up factors  
    dat_snp_all$accession <- factor(dat_snp_all$accession)
    dat_snp_all$section_block <- factor(dat_snp_all$section_block)
    dat_snp_all$section <- factor(dat_snp_all$section)
    dat_snp_all$locality <- factor(dat_snp_all$locality)
    
  # Scale PCA axes
    dat_snp_all$PC1_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC1_gen)) / sd(dat_snp_all$PC1_gen)
    dat_snp_all$PC2_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC2_gen)) / sd(dat_snp_all$PC2_gen)
    dat_snp_all$PC3_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC3_gen)) / sd(dat_snp_all$PC3_gen)
    dat_snp_all$PC4_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC4_gen)) / sd(dat_snp_all$PC4_gen)
    dat_snp_all$PC5_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC5_gen)) / sd(dat_snp_all$PC5_gen)
    
  # Save unscaled version  
    dat_snp_all_unscaled <- dat_snp_all
 
    summary(dat_snp_all[, 1050:1060])
    
  
  # Save means and sds for prediction
    scaled_snps_means_all <- apply(dat_snp_all[, snp_col_names], MARGIN = 2,
                                   function(x) mean(x, na.rm = TRUE))
    scaled_snps_sds_all <- apply(dat_snp_all[, snp_col_names], MARGIN = 2,
                                 function(x) sd(x, na.rm = TRUE))

    ## Scale genotypes
    dat_snp_all[, snp_col_names] <-   apply(dat_snp_all[, 
                                                snp_col_names], MARGIN = 2,
                                        function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))


 ## Create folds for cross validation
    
    
    # Select a subset of adults for taining and testing, CV Validation
    
    
    # 70/30 sampling
    # For path "run_475547_tmax_sum_dif_training_set_resids_7030_2019_02_26"
    # Used in initial submission to PNAS
    # set.seed(129)
    # training_moms <- sample(na.omit(climate_gbs_mom$accession),
    #                         size = length(na.omit(climate_gbs_mom$accession))*0.70)
    # 
    # testing_moms <- na.omit(climate_gbs_mom$accession)[!na.omit(climate_gbs_mom$accession)
    #                                                    %in% training_moms]
    # 
    # any(training_moms %in% testing_moms) # Should both be false
    # any(testing_moms %in% training_moms)
    # 
    # 
    # ## Save two datasets for training and testing
    # dat_gbs_training_clim <- dat_gbs_all_clim %>%
    #   dplyr::filter(accession %in% training_moms)
    # dat_gbs_testing_clim <- dat_gbs_all_clim %>%
    #   dplyr::filter(accession %in% testing_moms)
    # 
    # dim(dat_gbs_training_clim)
    # dim(dat_gbs_testing_clim)
    # 
    # table(dat_gbs_testing_clim$section_block)
    # table(dat_gbs_testing_clim$locality)
    # table(dat_gbs_testing_clim$accession)
    # 
    # table(dat_gbs_training_clim$section_block)
    # table(dat_gbs_training_clim$locality)
    # table(dat_gbs_training_clim$accession)
    # 
    # unique(dat_gbs_testing_clim$accession) %in% unique(dat_gbs_training_clim$accession)
    # unique(dat_gbs_training_clim$accession) %in% unique(dat_gbs_testing_clim$accession)
    # sum(unique(dat_gbs_training_clim$accession) %in% unique(dat_gbs_testing_clim$accession))
    
    
    # 10 fold cross validation - not across families
    set.seed(129)
    
    training_folds_moms <-  caret::createFolds(na.omit(climate_gbs_mom$accession),
                                          k = 10, returnTrain = T)
    
    training_folds <- list()
    for(fold in 1:length(training_folds_moms)){
      training_folds[[fold]] <- which(dat_gbs_all_clim$accession %in% 
                                  na.omit(climate_gbs_mom$accession)[training_folds_moms[[fold]]]) 
    }
    
    
    
    ## Check on folds
    x = 1
    for(fold in 1:length(training_folds)){
      cat(paste0("|--- Fold ", fold, " ---| \n"))
      cat("percent in training vs. testing: \n ")
      print(nrow(dat_gbs_all_clim[training_folds[[fold]], ])/nrow(dat_gbs_all_clim))
      print(nrow(dat_gbs_all_clim[-training_folds[[fold]], ])/nrow(dat_gbs_all_clim))
      
      cat("Number of accessions in training that are in testing \n ")
      print(sum(unique(dat_gbs_all_clim$accession[training_folds[[fold]]]) %in% 
                  unique(dat_gbs_all_clim$accession[-training_folds[[fold]]])))
      x = x + 1
    }

    
## 10 fold validation across families    
    # set.seed(129)
    # 
    # training_folds <-  caret::createFolds(factor(dat_gbs_all_clim$accession),
    #                                           k = 10, returnTrain = T)
    # 
    # ## Check on folds
    # x = 1
    # for(fold in 1:length(training_folds)){
    #   cat(paste0("|--- Fold ", fold, " ---| \n"))
    #   cat("percent in training vs. testing: \n ")
    #   print(nrow(dat_gbs_all_clim[training_folds[[fold]], ])/nrow(dat_gbs_all_clim))
    #   print(nrow(dat_gbs_all_clim[-training_folds[[fold]], ])/nrow(dat_gbs_all_clim))
    #   
    #   cat("Number of accessions in training that are in testing \n ")
    #   print(sum(unique(dat_gbs_all_clim$accession[training_folds[[fold]]]) %in% 
    #         unique(dat_gbs_all_clim$accession[-training_folds[[fold]]])))
    #   x = x + 1
    # }
    
  ### Run full model and save out residuals for cluster analysis
  
  ## Original used in PNAS submission
  fixed_effects <- paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"cr\") + s(tmax_sum_dif, bs=\"cr\") + s(accession, bs = \"re\") + s(locality, bs = \"re\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\")  + s(tmax_sum_dif, by = PC1_gen, bs=\"cr\") + s(tmax_sum_dif, by = PC2_gen, bs=\"cr\") + s(PC1_clim, bs=\"cr\") + s(PC2_clim, bs=\"cr\")  + s(tmax_sum_dif, by = PC1_clim, bs=\"cr\") + s(tmax_sum_dif, by = PC2_clim, bs=\"cr\")"))

  
  # Gam with interaction effect
  gam_snp_all = bam(formula = formula(fixed_effects),
                    data = dat_snp_all,
                    discrete = FALSE, 
                    nthreads = 8,
                    method = "fREML",
                    family = "tw",
                    control = list(trace = FALSE))
  
  summary(gam_snp_all)
  
  dat_snp_all$rgr_resids <- resid(gam_snp_all)
  dat_snp_all_unscaled$rgr_resids <- resid(gam_snp_all)
  hist(dat_snp_all$rgr_resids, breaks = 50)
  
  save.image(paste0("./output/gam_cluster_", Sys.Date(), ".Rdata"))
  
  
  summary(back_transform(dat_snp_all$tmax_sum_dif,
                         "tmax_sum_dif",
                         scaled_var_means_gbs_all,
                         scaled_var_sds_gbs_all))


 