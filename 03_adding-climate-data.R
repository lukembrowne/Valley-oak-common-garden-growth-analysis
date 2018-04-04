
## Requires 02_processing-climate-rasters.R to be run and files output


# Load libraries ----------------------------------------------------------

  library(psych)


# Load in climate data ----------------------------------------------------


  ## Read in climate data of of maternal trees - "Valley oak maternal tree climate data BCM 2018_03_08"
  ## If token is stale: gs_auth(new_user = TRUE)
 # climate_mom <- gs_read(gs_key("1lOQcBjdjUt-I4brnrD5LTCODZ2lvgzpqNnpxkYt_oSE"), ws = 2)
  
  ## Or read in 1921-1950 data
  climate_mom <- read_csv("./data/cleaned_data/GBS tree climate data BCM 1921-1950 2018_03_28.csv")

  ## Rename columns
  climate_mom <- climate_mom %>%
    dplyr::select(-Locality) %>%
    rename(Elevation = `Elevation (m)`) %>%
    select_all(., tolower)
  
  climate_mom
    
  ## Merge data
  dat_all_clim <- left_join(dat_all, climate_mom, by = "accession")
  dim(dat_all_clim)
  glimpse(dat_all_clim)
  
  ## Load in climate of common garden locations - 2010 to 2015
  garden_climate <- read_csv("./data/cleaned_data/common garden climate data BCM 2014-2016 2018_03_08.csv")
  
  glimpse(garden_climate)




# Calculating difference in climate variables -----------------------------

  ## Core variables from Riordan et al. 2016 Am J Botany  
  ## Except we are excluding AET because it is very strongly correlated with cwd
  climate_vars <- c("tmax_sum", "tmin_winter", "cwd", "bioclim_04", "bioclim_15",
                    "bioclim_18", "bioclim_19")
  
  climate_vars_dif <- paste(climate_vars, "_dif", sep = "")
  
  x = 1
  
  ## Calculate difference in climate as ## GARDEN - PROVENANCE
  for(var in climate_vars){
    
    dat_all_clim[dat_all_clim$site == "Chico", climate_vars_dif[x]] <- as.numeric(garden_climate[garden_climate$site == "Chico", var]) - dat_all_clim[dat_all_clim$site == "Chico", var]
    
    dat_all_clim[dat_all_clim$site == "IFG", climate_vars_dif[x]] <- as.numeric(garden_climate[garden_climate$site == "IFG", var]) - dat_all_clim[dat_all$site == "IFG", var]
    
    x = x + 1
  }
  
  
  hist(dat_all_clim$tmax_sum_dif)
  
  summary(dat_all_clim[, climate_vars_dif])
  
  ## Find NAs if there are any
  #View(dat_all_clim[is.na(pull(dat_all_clim, climate_vars_dif[1])), ])
  

  
  
## Correlations between climate variables
    
  ## Climate variables in seedlings
 # pairs.panels(dat_all_clim[, climate_vars], scale = TRUE)
  
  ## Climate variables for moms
 # pairs.panels(climate_mom[, climate_vars], scale = TRUE)
  
  ## Climate distances for seedlings
 # pairs.panels(dat_all_clim[, climate_vars_dif], scale = TRUE)



# Scaling predictor variables -------------------------------------------------------

  dat_all_scaled <- dat_all_clim
  
  x = 1
  
  scaled_var_means <- NA
  scaled_var_sds <- NA
  
  to_scale = colnames(dplyr::select(dat_all_clim, climate_vars_dif, 
                                    height_2015, height_2014,
                                    climate_vars))
  
  for(var in to_scale){
    
    scaled_var_means[x] <- mean(dat_all_scaled[[var]], na.rm = TRUE)
    scaled_var_sds[x] <- sd(dat_all_scaled[[var]], na.rm = TRUE)
    
    dat_all_scaled[var] <- (dat_all_scaled[var] - scaled_var_means[x]) / scaled_var_sds[x]
    
    x = x + 1
  }
  
  names(scaled_var_means) <- to_scale
  names(scaled_var_sds) <- to_scale
  
  
  summary(dat_all_scaled[, to_scale]) ## Means should all be 0
  



# 
# # PCA on climate variables ------------------------------------------------
# 
#   
#   
#   library(factoextra)
#   
#   clim_pca <- prcomp(dplyr::select(dat_all_scaled, climate_vars_dif), scale = FALSE) 
#   
#   summary(clim_pca)
# 
#   
#   ## PCA on home variables
#   clim_pca_home <- prcomp(dat_all_scaled[, climate_vars], scale = FALSE) 
#   
#   summary(clim_pca_home)
#   
#   
#   ## Format as dataframe
#   clim_pca_vals <- as.data.frame(clim_pca$x)
#   colnames(clim_pca_vals) <- paste(colnames(clim_pca_vals), "_clim_dif", sep = "")
#   
#   clim_pca_home_vals <- as.data.frame(clim_pca_home$x)
#   colnames(clim_pca_home_vals) <- paste(colnames(clim_pca_home_vals), "_home", sep = "")
#   
#   
#   ## Scree plot of PCA
#   fviz_screeplot(clim_pca)
#   
#   
#   ## Biplot of PCA
#   # autoplot(clim_pca, data = dat_all_scaled, colour = "prov", loadings = TRUE, 
#   #          loadings.label = TRUE) + theme_bw() + theme(legend.position="none") 
#   
#   
#   ## Contributions of each axis
#   # fviz_contrib(clim_pca, choice = "var", axes = c(1))
#   # fviz_contrib(clim_pca, choice = "var", axes = c(2))
#   # fviz_contrib(clim_pca, choice = "var", axes = c(3))
#   # fviz_contrib(clim_pca, choice = "var", axes = c(4))
#   # fviz_contrib(clim_pca, choice = "var", axes = c(5))
#   # 
#   ## Arrows plot
#   fviz_pca_var(clim_pca, axes = c(1,2))
#   
#   
#   ## Join back to main dataframe
#   dat_all_scaled <- bind_cols(dat_all_scaled, clim_pca_vals)


