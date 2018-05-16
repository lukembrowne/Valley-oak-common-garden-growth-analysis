
## Requires 02_processing-climate-rasters.R to be run and files output


# Load libraries ----------------------------------------------------------

  library(psych)


# Load in climate data ----------------------------------------------------

  ## Read in climate data of of maternal trees in the common garden
  ## Should be 659 trees
  climate_garden_mom <- read_csv("./data/cleaned_data/maternal tree climate data BCM 1950-1981 2018_03_08.csv")
  
  ## Rename columns
  climate_garden_mom <- climate_garden_mom %>%
    dplyr::select(-Locality) %>%
    rename(Elevation = `Elevation (m)`) %>%
    select_all(., tolower)
  
  climate_garden_mom
  
  ## Read in climate data for GBS trees
  climate_gbs_mom <- read_csv("./data/cleaned_data/GBS tree climate data BCM 1950-1981 2018_03_08.csv")
  
  climate_gbs_mom <- climate_gbs_mom %>%
    select_all(., tolower)
  
  climate_gbs_mom 
  
  climate_gbs_mom <- left_join(climate_gbs_mom, 
                               dplyr::select(gen_dat, gbs_name, accession),
                               by = c("id" = "gbs_name"))
  
  climate_gbs_mom$accession
  
  ## Merge data to give climate of origin for each seedling
    dat_all_clim <- left_join(dat_all, climate_garden_mom, by = "accession")
    glimpse(dat_all_clim)
    dim(dat_all_clim)
  
    
  ## Load in climate of common garden locations - 2010 to 2015
    garden_climate <- read_csv("./data/cleaned_data/common garden climate data BCM 2014-2016 2018_03_08.csv")
    glimpse(garden_climate)

    
# Filter down to seedlings with and without GBS data ----------------------
    
    ## Filtering down individuals who's moms don't have climate data and has GBS data
    dat_gbs_only_clim <- dplyr::filter(dat_all_clim, accession %in% climate_gbs_mom$accession)
    dim(dat_gbs_only_clim)


# Choose climate variables ------------------------------------------------
   
     ## Core variables from Riordan et al. 2016 Am J Botany  
    ## Except we are excluding AET because it is very strongly correlated with cwd
    climate_vars <- c("tmax_sum",
                      "tmin_winter",
                      #  "tmax",
                      #"tmin",
                      # "latitude",
                      #  "longitude",
                      #  "elevation")
                      #"cwd", 
                      "bioclim_04") # Temperature seasonality
    #"bioclim_15", # Precipitation seasonality
    # "bioclim_18", # Precipitation of Warmest Quarter
    #"bioclim_19") # Precipitation of Coldest Quarter
    
    
    
# PCA on climate variables ------------------------------------------------

    library(factoextra)

    clim_pca_home <- prcomp(dat_all_clim[, climate_vars], scale = TRUE)
    summary(clim_pca_home)
    
    
    ## Format as dataframe
    clim_pca_home_vals <- as.data.frame(clim_pca_home$x)
    
    ## Scree plot of PCA
    fviz_screeplot(clim_pca_home)
    
    ## Contributions of each axis
    fviz_contrib(clim_pca_home, choice = "var", axes = c(1))
    fviz_contrib(clim_pca_home, choice = "var", axes = c(2))
    fviz_contrib(clim_pca_home, choice = "var", axes = c(3))
    fviz_contrib(clim_pca_home, choice = "var", axes = c(4))
    fviz_contrib(clim_pca_home, choice = "var", axes = c(5))
    
    ## Arrows plot
    fviz_pca_var(clim_pca_home, axes = c(1,2))
    
    ## Join back to main dataframe
    dat_all_clim <- bind_cols(dat_all_clim, clim_pca_home_vals)
    
    ## Predict PC values at Chico  & IFG based on pca of climate of origin
    # Need to be sure that the garden climate variables are scaled first
    
    garden_climate_scaled <- garden_climate # Scaled based on seedling climate of origin!!!
    
    for(var in climate_vars){
      garden_climate_scaled[, var] <- (pull(garden_climate_scaled, var) - mean(pull(dat_all_clim, var))) / sd(pull(dat_all_clim, var))
    }
    
    garden_climate_scaled[, climate_vars]
    
    garden_climate_pca_vals <-   as.data.frame(predict(clim_pca_home, garden_climate_scaled))
    garden_climate_pca_vals
    
    garden_climate <- bind_cols(garden_climate, garden_climate_pca_vals)
    
    ## Add PCs to climate vars
      climate_vars <- c(climate_vars, colnames(garden_climate_pca_vals))
    
    
# Calculating difference in climate variables -----------------------------

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
  

  

# Correlations in climate variables ---------------------------------------

  ## Climate variables in seedlings
  pairs.panels(dat_all_clim[, climate_vars], scale = TRUE)

  ## Climate variables for moms
 # pairs.panels(climate_mom[, climate_vars], scale = TRUE)
  
  ## Climate distances for seedlings
  pairs.panels(dat_all_clim[, climate_vars_dif], scale = TRUE)


# Scaling predictor variables -------------------------------------------------------
  
  # Need to run separately for both full dataset and dataset with only gbs samples

  dat_all_scaled <- dat_all_clim
  dat_gbs_only_scaled <- dat_gbs_only_clim
  
  x = 1
  
  scaled_var_means_gbs_only <- NA
  scaled_var_sds_gbs_only <- NA
  
  scaled_var_means_all <- NA
  scaled_var_sds_all <- NA
  
  to_scale = colnames(dplyr::select(dat_all_clim, climate_vars_dif, 
                                    height_2015, height_2014,
                                    climate_vars))
  
  for(var in to_scale){
    
    # For all seedlings  
      scaled_var_means_all[x] <- mean(dat_all_scaled[[var]], na.rm = TRUE)
      scaled_var_sds_all[x] <- sd(dat_all_scaled[[var]], na.rm = TRUE)
      dat_all_scaled[var] <- (dat_all_scaled[var] - scaled_var_means_all[x]) / scaled_var_sds_all[x]
      
    # For seedlings with gbs data  
      scaled_var_means_gbs_only[x] <- mean(dat_gbs_only_scaled[[var]], na.rm = TRUE)
      scaled_var_sds_gbs_only[x] <- sd(dat_gbs_only_scaled[[var]], na.rm = TRUE)
      dat_gbs_only_scaled[var] <- (dat_gbs_only_scaled[var] - scaled_var_means_gbs_only[x]) /
                                  scaled_var_sds_gbs_only[x]
  
      x = x + 1
  }
  
# Assign names  
  names(scaled_var_means_all) <- to_scale
  names(scaled_var_sds_all) <- to_scale
  
  names(scaled_var_means_gbs_only) <- to_scale
  names(scaled_var_sds_gbs_only) <- to_scale
  
# Check out results
  scaled_var_means_all
  scaled_var_means_gbs_only
  
  scaled_var_sds_all
  scaled_var_sds_gbs_only
  
  summary(dat_all_scaled[, to_scale]) ## Means should all be 0
  
  summary(dat_gbs_only_scaled[, to_scale]) ## Means should all be 0
  




  
  
  
  
  
  
  
    


