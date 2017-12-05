
# Scaling climate data prior to PCA - grouping by mom since all progeny will have the same climate

dat_all_scaled_mom <- dat_all %>%
  group_by(mom) %>%
  summarise_all(mean, na.rm = TRUE)

x = 1

climate_var_means <- NA
climate_var_sds <- NA

for(var in c(climate_vars_dif, climate_vars)){
  
  climate_var_means[x] <- mean(dat_all_scaled_mom[[var]])
  climate_var_sds[x] <- sd(dat_all_scaled_mom[[var]])
  
  dat_all_scaled_mom[var] <- (dat_all_scaled_mom[var] - climate_var_means[x]) / climate_var_sds[x]
  
  x = x+1
}

names(climate_var_means) <- c(climate_vars_dif, climate_vars)
names(climate_var_sds) <- c(climate_vars_dif, climate_vars)


## Run PCA on climate variables
clim_pca <- prcomp(dplyr::select(dat_all_scaled_mom, contains("_dif")), scale = FALSE) 

summary(clim_pca)


## PCA on home variables
clim_pca_home <- prcomp(dat_all_scaled_mom[, climate_vars], scale = FALSE) 

summary(clim_pca_home)



## Format as dataframe
clim_pca_vals <- as.data.frame(clim_pca$x)
clim_pca_home_vals <- as.data.frame(clim_pca_home$x)
colnames(clim_pca_home_vals) <- paste(colnames(clim_pca_home_vals), "_home", sep = "")

## Scree plot of PCA
fviz_screeplot(clim_pca)

## Biplot of PCA
autoplot(clim_pca, data = dat_all_scaled_mom, colour = "prov", loadings = TRUE, 
         loadings.label = TRUE) + theme_bw() + theme(legend.position="none") 

## Contributions of each axis
fviz_contrib(clim_pca, choice = "var", axes = c(1))
fviz_contrib(clim_pca, choice = "var", axes = c(2))
fviz_contrib(clim_pca, choice = "var", axes = c(3))
fviz_contrib(clim_pca, choice = "var", axes = c(4))
fviz_contrib(clim_pca, choice = "var", axes = c(5))

## Arrows plot
fviz_pca_var(clim_pca, axes = c(1,2))

dat_all_pca <- left_join(dat_all, bind_cols(dat_all_scaled_mom[, "mom"], clim_pca_vals))
# dat_all_pca <- left_join(dat_all_pca, bind_cols(dat_all_scaled_mom[, "mom"],
#                                             clim_pca_home_vals))


### SHOULD THIS BE SCALED BASED ON MOMS??
elev_sd <- sd(dat_all_pca$elev)
elev_mean <- mean(dat_all_pca$elev)

dat_all_pca$elev <- (dat_all_pca$elev - elev_mean) / elev_sd

dat_all_pca$height_2016_scaled <- (dat_all_pca$height_2016 - mean(dat_all_pca$height_2016, na.rm = T)) / sd(dat_all_pca$height_2016, na.rm = T)

