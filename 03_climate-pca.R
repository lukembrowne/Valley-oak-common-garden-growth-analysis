

# Load libraries ----------------------------------------------------------

library(factoextra)


# PCA on climate variables ------------------------------------------------

  clim_pca <- prcomp(dplyr::select(dat_all_scaled, contains("_dif")), scale = FALSE) 
  
  summary(clim_pca)

  
  ## PCA on home variables
  clim_pca_home <- prcomp(dat_all_scaled[, climate_vars], scale = FALSE) 
  
  summary(clim_pca_home)
  
  
  ## Format as dataframe
  clim_pca_vals <- as.data.frame(clim_pca$x)
  colnames(clim_pca_vals) <- paste(colnames(clim_pca_vals), "_clim_dif", sep = "")
  
  clim_pca_home_vals <- as.data.frame(clim_pca_home$x)
  colnames(clim_pca_home_vals) <- paste(colnames(clim_pca_home_vals), "_home", sep = "")
  
  
  ## Scree plot of PCA
  fviz_screeplot(clim_pca)
  
  
  ## Biplot of PCA
  ggplot2::autoplot(clim_pca, data = dat_all_scaled, colour = "prov", loadings = TRUE, 
           loadings.label = TRUE) + theme_bw() + theme(legend.position="none") 
  
  
  ## Contributions of each axis
  # fviz_contrib(clim_pca, choice = "var", axes = c(1))
  # fviz_contrib(clim_pca, choice = "var", axes = c(2))
  # fviz_contrib(clim_pca, choice = "var", axes = c(3))
  # fviz_contrib(clim_pca, choice = "var", axes = c(4))
  # fviz_contrib(clim_pca, choice = "var", axes = c(5))
  # 
  ## Arrows plot
  fviz_pca_var(clim_pca, axes = c(1,2))
  
  
  ## Join back to main dataframe
  dat_all_scaled <- bind_cols(dat_all_scaled, clim_pca_vals)


