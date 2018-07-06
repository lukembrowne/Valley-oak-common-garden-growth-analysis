
## Requires 02_processing-climate-rasters.R to be run and files output


# Load libraries ----------------------------------------------------------

  library(psych)
  library(factoextra)

# Load in climate data ----------------------------------------------------

  ## Read in BCM climate data of of maternal trees in the common garden
  ## Should be 659 trees
    climate_garden_mom <- read_csv("./data/cleaned_data/maternal tree climate data BCM 1950-1981 2018_03_08.csv")
  
  ## Rename columns
    climate_garden_mom <- climate_garden_mom %>%
      dplyr::select(-Locality) %>%
      rename(Elevation = `Elevation (m)`) %>%
      select_all(., tolower)
    
    climate_garden_mom
    
 
  ## Function for calculating growing degree days     
  # Code for GDD - https://rdrr.io/cran/pollen/src/R/gdd.R  
    
    test_0 <- function(x){
      ifelse(x < 0, 0, x)
    }

   # tmax = climate_garden_mom$tmax_jan
    #tmin = climate_garden_mom$tmin_jan
    tbase = 5
    tbase_max = 50


    adjust_for_tbase <- function(x, tbase) {
      ifelse(test = x < tbase, yes = tbase, no = x)
    }
    adjust_for_tbase_max <- function(x, tbase_max) {
      ifelse(test = x > tbase_max, yes = tbase_max, no = x)
    }

    # tmax_adjusted <- adjust_for_tbase(tmax, tbase)
    # tmin_adjusted <- adjust_for_tbase(tmin, tbase)
    # 
    # tmax_adjusted <- adjust_for_tbase_max(tmax_adjusted, tbase_max)
    # tmin_adjusted <- adjust_for_tbase_max(tmin_adjusted, tbase_max)
    # 
    # gdd_temp <- (tmax_adjusted + tmin_adjusted) / 2 - tbase
    # 
  calc_gdd <- function(dat, threshold){
    
      (  test_0(((adjust_for_tbase_max(dat$tmax_jan, tbase_max) + dat$tmin_jan) / 2 - threshold)) * 31) + # Jan
      (  test_0(((adjust_for_tbase_max(dat$tmax_feb, tbase_max) + dat$tmin_feb) / 2 - threshold)) * 28) + # Feb
      (  test_0(((adjust_for_tbase_max(dat$tmax_mar, tbase_max) + dat$tmin_mar) / 2 - threshold)) * 31) + # Mar
      (  test_0(((adjust_for_tbase_max(dat$tmax_apr, tbase_max) + dat$tmin_apr) / 2 - threshold)) * 30) + # Apr
      (  test_0(((adjust_for_tbase_max(dat$tmax_may, tbase_max) + dat$tmin_may) / 2 - threshold)) * 31) + # May
      (  test_0(((adjust_for_tbase_max(dat$tmax_jun, tbase_max) + dat$tmin_jun) / 2 - threshold)) * 30) + # Jun
      (  test_0(((adjust_for_tbase_max(dat$tmax_jul, tbase_max) + dat$tmin_jul) / 2 - threshold)) * 31) + # Jul
      (  test_0(((adjust_for_tbase_max(dat$tmax_aug, tbase_max) + dat$tmin_aug) / 2 - threshold)) * 31) + # Aug
      (  test_0(((adjust_for_tbase_max(dat$tmax_sep, tbase_max) + dat$tmin_sep) / 2 - threshold)) * 30) + # Sep
      (  test_0(((adjust_for_tbase_max(dat$tmax_oct, tbase_max) + dat$tmin_oct) / 2 - threshold)) * 31) + # Oct
      (  test_0(((adjust_for_tbase_max(dat$tmax_nov, tbase_max) + dat$tmin_nov) / 2 - threshold)) * 30) + # Nov
      (  test_0(((adjust_for_tbase_max(dat$tmax_dec, tbase_max) + dat$tmin_dec) / 2 - threshold)) * 31)   # Dec
    }
      
  climate_garden_mom$DD5 <- calc_gdd(dat = climate_garden_mom, threshold = 5)
    
  summary(climate_garden_mom$DD5)
    
  ## Read in historical / paleo climate data from WNA to compare with BCM
    
    # Last 1000 years
      last1000 <-   read_csv("./data/gis/climate_data/climateWNA/maternal tree climate data ClimateWNA Last1000 2018_06_11_4GCM-Ensemble_past1000MSY.csv") %>%
                    rename(accession = ID2 , 
                           tmax_sum = Tmax_sm,
                           tmin_winter = Tmin_wt,
                           tave = MAT) %>%
                    dplyr::select(-ID1, -Latitude, -Longitude, -Elevation)
      
      colnames(last1000)[-1] <- paste0(colnames(last1000)[-1], "_last1000")
      
    # Last glacial maximum
      lgm <- read_csv("./data/gis/climate_data/climateWNA/maternal tree climate data ClimateWNA LGM 2018_06_11_4GCM-Ensemble_lgmMSY.csv") %>%
        rename(accession = ID2 , 
               tmax_sum = Tmax_sm,
               tmin_winter = Tmin_wt,
               tave = MAT) %>%
        dplyr::select(-ID1, -Latitude, -Longitude, -Elevation)
      
      colnames(lgm)[-1] <- paste0(colnames(lgm)[-1], "_lgm")
    
    # Mid-holocene  
      holo <- read_csv("./data/gis/climate_data/climateWNA/maternal tree climate data ClimateWNA Holocene 2018_06_11_4GCM-Ensemble_midHoloceneMSY.csv") %>%
        rename(accession = ID2 , 
               tmax_sum = Tmax_sm,
               tmin_winter = Tmin_wt,
               tave = MAT) %>%
        dplyr::select(-ID1, -Latitude, -Longitude, -Elevation)
      
      colnames(holo)[-1] <- paste0(colnames(holo)[-1], "_holo")
      
      
    # Join together in main climate dataset
      climate_garden_mom <- left_join(left_join(left_join(climate_garden_mom, last1000), lgm), holo)
      
      dim(climate_garden_mom)
      
      
  ## Read in climate data for GBS trees
    climate_gbs_mom <- read_csv("./data/cleaned_data/GBS tree climate data BCM 1950-1981 2018_03_08.csv")
    
    climate_gbs_mom <- climate_gbs_mom %>%
      select_all(., tolower)
    
    climate_gbs_mom 
    
    climate_gbs_mom <- left_join(climate_gbs_mom, 
                                 dplyr::select(gen_dat, gbs_name, accession),
                                 by = c("id" = "gbs_name"))
    
    climate_gbs_mom$accession
    
    # Reading in paleo climate data
    # Last 1000 years
    last1000 <-   read_csv("./data/gis/climate_data/climateWNA/GBS tree climate data Climate WNA 2018_06_12_4GCM-Ensemble_past1000MSY.csv") %>%
      rename(id = ID1 , 
             tmax_sum = Tmax_sm,
             tmin_winter = Tmin_wt,
             tave = MAT) %>%
      dplyr::select(-ID2, -Latitude, -Longitude, -Elevation)
    
    colnames(last1000)[-1] <- paste0(colnames(last1000)[-1], "_last1000")
    
    # Last glacial maximum
    lgm <- read_csv("./data/gis/climate_data/climateWNA/GBS tree climate data Climate WNA 2018_06_12_4GCM-Ensemble_lgmMSY.csv") %>%
      rename(id = ID1 , 
             tmax_sum = Tmax_sm,
             tmin_winter = Tmin_wt,
             tave = MAT) %>%
      dplyr::select(-ID2, -Latitude, -Longitude, -Elevation)
    
    colnames(lgm)[-1] <- paste0(colnames(lgm)[-1], "_lgm")
    
    # Mid-holocene  
    holo <- read_csv("./data/gis/climate_data/climateWNA/GBS tree climate data Climate WNA 2018_06_12_4GCM-Ensemble_midHoloceneMSY.csv") %>%
      rename(id = ID1 , 
             tmax_sum = Tmax_sm,
             tmin_winter = Tmin_wt,
             tave = MAT) %>%
      dplyr::select(-ID2, -Latitude, -Longitude, -Elevation)
    colnames(holo)[-1] <- paste0(colnames(holo)[-1], "_holo")
    
    
    # Join together in main climate dataset
    climate_gbs_mom <- left_join(left_join(left_join(climate_gbs_mom, last1000), lgm), holo)
    
    dim(climate_gbs_mom)
      
  
      
  ## Add random climate variable for testing
    climate_gbs_mom$random <- rnorm(n = nrow(climate_gbs_mom))
    climate_garden_mom$random <- rnorm(n = nrow(climate_garden_mom))
    
    
    ## Calculate average temperature
    climate_gbs_mom$tave <- (climate_gbs_mom$tmax + climate_gbs_mom$tmin) / 2 
    climate_garden_mom$tave <- (climate_garden_mom$tmax + climate_garden_mom$tmin) / 2 
    
    
    
  ## Calculating MEMs
    
    library(spmoran)
    
    # Moms in common garden
    mem <- meigen(climate_garden_mom[, c("longitude", "latitude")])
    
    mem$ev # Eigenvalues
    mem$sf # Eigenvectors
    
    climate_garden_mom$mem1 <- mem$sf[, 1]
    climate_garden_mom$mem2 <- mem$sf[, 2]
    climate_garden_mom$mem3 <- mem$sf[, 3]
    climate_garden_mom$mem4 <- mem$sf[, 4]
    
    # ggplot(climate_garden_mom, aes(x = longitude, y = latitude, size = mem1, col = mem1)) +
    # geom_point() 
    
    
    # GBS samples
    mem <- meigen(climate_gbs_mom[, c("longitude", "latitude")])
    
    mem$ev # Eigenvalues
    mem$sf # Eigenvectors
    
    climate_gbs_mom$mem1 <- mem$sf[, 1]
    climate_gbs_mom$mem2 <- mem$sf[, 2]
    climate_gbs_mom$mem3 <- mem$sf[, 3]
    climate_gbs_mom$mem4 <- mem$sf[, 4]
    
    
  ## Merge data to give climate of origin for each seedling
    dat_all_clim <- left_join(dat_all, climate_garden_mom, by = "accession")
    glimpse(dat_all_clim)
    dim(dat_all_clim)
  
    
  ## Load in climate of common garden locations - 2014-2016
    garden_climate <- read_csv("./data/cleaned_data/common garden climate data BCM 2014-2016 2018_03_08.csv")
    glimpse(garden_climate)
    
    # Add random variable
    garden_climate$random <- rnorm(n = nrow(garden_climate))
    
    # Add average temps
    garden_climate$tave <- (garden_climate$tmax + garden_climate$tmin) / 2
    
    # Add comparable variables for historic WNA
    garden_climate$tmax_sum_last1000 <- garden_climate$tmax_sum
    garden_climate$tmax_sum_lgm <- garden_climate$tmax_sum
    garden_climate$tmax_sum_holo <- garden_climate$tmax_sum
    
    garden_climate$tmin_winter_last1000 <- garden_climate$tmin_winter
    garden_climate$tmin_winter_lgm <- garden_climate$tmin_winter
    garden_climate$tmin_winter_holo <- garden_climate$tmin_winter
    
    garden_climate$tave_last1000 <- garden_climate$tave
    garden_climate$tave_lgm <- garden_climate$tave
    garden_climate$tave_holo <- garden_climate$tave
    
    ## Calculating growing degree days
    garden_climate$DD5 <- calc_gdd(dat = garden_climate, threshold = 5)
    garden_climate$DD5_lgm <- garden_climate$DD5
    
    
    
# Choose climate variables ------------------------------------------------
   
     ## Core variables from Riordan et al. 2016 Am J Botany  
    ## Except we are excluding AET because it is very strongly correlated with cwd
    climate_vars <- c("tmax_sum",
                      "tmax","tmin",
                      "tave",
                      "tave_lgm",
                    #  "tmax_sum_last1000",
                      "tmax_sum_lgm",
                    #  "tmax_sum_holo",
                      "tmin_winter",
                    #  "tmin_winter_last1000",
                      "tmin_winter_lgm",
                     # "tmin_winter_holo",
                      "DD5",
                    "DD5_lgm",
                  #  "bioclim_04",
                      "random")
                      #  "tmax",
                      #"tmin",
                      # "latitude",
                      #  "longitude",
                      #  "elevation")
                      #"cwd", 
                      #"bioclim_04") # Temperature seasonality
    #"bioclim_15", # Precipitation seasonality
    # "bioclim_18", # Precipitation of Warmest Quarter
    #"bioclim_19") # Precipitation of Coldest Quarter
    
    
    
# PCA on climate variables ------------------------------------------------

  #   # PCA based on climate data of maternal trees in provenance trial
  #     clim_pca_home <- prcomp(climate_garden_mom[,
  #                                                climate_vars[-which(climate_vars == "random")]], scale = TRUE)
  #     summary(clim_pca_home)
  # 
  # 
  #   ## Format as dataframe
  #     clim_pca_home_vals <- as.data.frame(clim_pca_home$x)
  # 
  #   ## Scree plot of PCA
  #     fviz_screeplot(clim_pca_home)
  # 
  #   ## Contributions of each axis
  #     # fviz_contrib(clim_pca_home, choice = "var", axes = c(1))
  #     # fviz_contrib(clim_pca_home, choice = "var", axes = c(2))
  #     # fviz_contrib(clim_pca_home, choice = "var", axes = c(3))
  #     # fviz_contrib(clim_pca_home, choice = "var", axes = c(4))
  #     # fviz_contrib(clim_pca_home, choice = "var", axes = c(5))
  # 
  #   ## Arrows plot
  #     fviz_pca_var(clim_pca_home, axes = c(1,2))
  # 
  #   ## Predict PC values for
  #   # 1) Seedlings in common garden
  #   # 2) Chico  & IFG based on pca of climate of origin
  #   # 3) GBS moms
  #   # Need to be sure that the garden climate variables are scaled first
  # 
  #   dat_all_clim_scaled <- dat_all_clim
  #   garden_climate_scaled <- garden_climate # Scaled based on seedling climate of origin!!!
  #   climate_gbs_mom_scaled <- climate_gbs_mom
  # 
  #   for(var in climate_vars){
  # 
  #     dat_all_clim_scaled[, var] <- (pull(dat_all_clim_scaled, var) - mean(pull(dat_all_clim, var))) / sd(pull(dat_all_clim, var))
  # 
  #     garden_climate_scaled[, var] <- (pull(garden_climate_scaled, var) - mean(pull(dat_all_clim, var))) / sd(pull(dat_all_clim, var))
  # 
  #    climate_gbs_mom_scaled[, var] <- (pull( climate_gbs_mom_scaled, var) - mean(pull(dat_all_clim, var))) / sd(pull(dat_all_clim, var))
  # 
  #   }
  # 
  #   dat_all_clim_scaled[, climate_vars]
  #   garden_climate_scaled[, climate_vars]
  #   climate_gbs_mom_scaled[, climate_vars]
  # 
  # # Predict PCA values
  # 
  #   dat_all_clim_pca_vals <-   as.data.frame(predict(clim_pca_home, dat_all_clim_scaled))
  #   dat_all_clim_pca_vals
  # 
  #   garden_climate_pca_vals <-   as.data.frame(predict(clim_pca_home, garden_climate_scaled))
  #   garden_climate_pca_vals
  # 
  #   climate_gbs_mom_pca_vals <- as.data.frame(predict(clim_pca_home, climate_gbs_mom_scaled))
  #   climate_gbs_mom_pca_vals
  # 
  # # Join back to main dataframes
  #   dat_all_clim <- bind_cols(dat_all_clim, dat_all_clim_pca_vals)
  #   garden_climate <- bind_cols(garden_climate, garden_climate_pca_vals)
  #   # climate_gbs_mom <- bind_cols(climate_gbs_mom, climate_gbs_mom_pca_vals)
  # 
  # 
  # 
  # ## Add PCs to list of climate vars
  #    climate_vars <- c(climate_vars, colnames(garden_climate_pca_vals))
  # 

      
# Filter down to seedlings with and without GBS data ----------------------
      
      ## Filtering down individuals who's moms don't have climate data and has GBS data
      dat_gbs_only_clim <- dplyr::filter(dat_all_clim,
                                         accession %in% climate_gbs_mom$accession)
      dim(dat_gbs_only_clim)
      
    
# Calculating difference in climate variables -----------------------------

  climate_vars_dif <- paste(climate_vars, "_dif", sep = "")
  
  x = 1
  
  ## Calculate difference in climate as ## GARDEN - PROVENANCE
    for(var in climate_vars){
      
      
    # For all seedlings    
      dat_all_clim[dat_all_clim$site == "Chico", climate_vars_dif[x]] <- as.numeric(garden_climate[garden_climate$site == "Chico", var]) - dat_all_clim[dat_all_clim$site == "Chico", var]
      
      dat_all_clim[dat_all_clim$site == "IFG", climate_vars_dif[x]] <- as.numeric(garden_climate[garden_climate$site == "IFG", var]) - dat_all_clim[dat_all_clim$site == "IFG", var]
      
      
    # For GBS seedlings  
      dat_gbs_only_clim[dat_gbs_only_clim$site == "Chico", climate_vars_dif[x]] <- as.numeric(garden_climate[garden_climate$site == "Chico", var]) - dat_gbs_only_clim[dat_gbs_only_clim$site == "Chico", var]
      
      dat_gbs_only_clim[dat_gbs_only_clim$site == "IFG", climate_vars_dif[x]] <- as.numeric(garden_climate[garden_climate$site == "IFG", var]) - dat_gbs_only_clim[dat_gbs_only_clim$site == "IFG", var]
      
      x = x + 1
    }
  
  
  hist(dat_all_clim$tmax_sum_dif)
  
  summary(dat_all_clim[, climate_vars_dif])
  
  ## Find NAs if there are any
  #View(dat_all_clim[is.na(pull(dat_all_clim, climate_vars_dif[1])), ])
  


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
  


 pairs.panels(climate_garden_mom[, c("latitude", "longitude", 
                                     #"mem1", "mem2", "mem3", 
                                     climate_vars)])
 
 pairs.panels(climate_garden_mom[, c("tmax_sum", "tmin_winter", "DD5", "tmax_sum_lgm","tmin_winter_lgm", "DD5_lgm")], ellipses = F)
 
 pairs.panels(climate_garden_mom[, c("tmax_sum",
                                     "tmax_sum_lgm",
                                     "tmin_winter", 
                                     "tmin_winter_lgm",
                                     "DD5", 
                                     "DD5_lgm")], ellipses = F, cex = 2)
 
 pairs.panels(climate_garden_mom[, c("tmax_sum",
                                     "tmax",
                                     "DD5",
                                     "tave",
                                     "bioclim_04")])
 
 #pairs.panels(dat_all_scaled[, climate_vars_dif])
 



