
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
    
  ## Read in current / historical / paleo climate data from WNA to compare with BCM
    
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
    
    # Reading in ClimateWNA paleo climate data
    
    # 1951-1980
    wna_current <-   read_csv("./data/gis/climate_data/climateWNA/GBS tree climate data Climate WNA 2018_06_12_Normal_1951_1980MSY.csv") %>%
      rename(id = ID1 , 
             tmax_sum = Tmax_sm,
             tmin_winter = Tmin_wt,
             tave = MAT,
             cwd = CMD) %>%
      dplyr::select(-ID2, -Latitude, -Longitude, -Elevation)
    
    # Calculate tmax and tmin
    wna_current$tmax <- rowMeans(dplyr::select(wna_current, Tmax01,Tmax02,Tmax03,Tmax04,Tmax05,Tmax06,Tmax07,Tmax08,Tmax09,Tmax10,Tmax11,Tmax12))
    wna_current$tmin <- rowMeans(dplyr::select(wna_current, Tmin01,Tmin02,Tmin03,Tmin04,Tmin05,Tmin06,Tmin07,Tmin08,Tmin09,Tmin10,Tmin11,Tmin12))
    
    colnames(wna_current)[-1] <- paste0(colnames(wna_current)[-1], "_wna_current")
    
    # Calculate bioclim variables
    bioclim_wna_current = climates::bioclim2(tmin = as.data.frame(dplyr::select(wna_current, Tmin01_wna_current,
                                                                        Tmin02_wna_current,
                                                                        Tmin03_wna_current,
                                                                        Tmin04_wna_current,
                                                                        Tmin05_wna_current,
                                                                        Tmin06_wna_current,
                                                                        Tmin07_wna_current,
                                                                        Tmin08_wna_current,
                                                                        Tmin09_wna_current,
                                                                        Tmin10_wna_current,
                                                                        Tmin11_wna_current,
                                                                        Tmin12_wna_current)),
                                     ## Remove site column
                                     tmax = as.data.frame(dplyr::select(wna_current, Tmax01_wna_current,
                                                                        Tmax02_wna_current,
                                                                        Tmax03_wna_current,
                                                                        Tmax04_wna_current,
                                                                        Tmax05_wna_current,
                                                                        Tmax06_wna_current,
                                                                        Tmax07_wna_current,
                                                                        Tmax08_wna_current,
                                                                        Tmax09_wna_current,
                                                                        Tmax10_wna_current,
                                                                        Tmax11_wna_current,
                                                                        Tmax12_wna_current)),
                                     prec = as.data.frame(dplyr::select(wna_current, PPT01_wna_current,
                                                                        PPT02_wna_current,
                                                                        PPT03_wna_current,
                                                                        PPT04_wna_current,
                                                                        PPT05_wna_current,
                                                                        PPT06_wna_current,
                                                                        PPT07_wna_current,
                                                                        PPT08_wna_current,
                                                                        PPT09_wna_current,
                                                                        PPT10_wna_current,
                                                                        PPT11_wna_current,
                                                                        PPT12_wna_current)),
                                     files.as.inputs = FALSE)
    
    
    colnames(bioclim_wna_current) <- paste0(colnames(bioclim_wna_current), "_", "wna_current")
    wna_current <- cbind(wna_current, bioclim_wna_current)
    
    
    # Last 1000 years
    last1000 <-   read_csv("./data/gis/climate_data/climateWNA/GBS tree climate data Climate WNA 2018_06_12_4GCM-Ensemble_past1000MSY.csv") %>%
      rename(id = ID1 , 
             tmax_sum = Tmax_sm,
             tmin_winter = Tmin_wt,
             tave = MAT,
             cwd = CMD) %>%
      dplyr::select(-ID2, -Latitude, -Longitude, -Elevation)
    
    colnames(last1000)[-1] <- paste0(colnames(last1000)[-1], "_last1000")
    
    # Last glacial maximum
    lgm <- read_csv("./data/gis/climate_data/climateWNA/GBS tree climate data Climate WNA 2018_06_12_4GCM-Ensemble_lgmMSY.csv") %>%
      rename(id = ID1 , 
             tmax_sum = Tmax_sm,
             tmin_winter = Tmin_wt,
             tave = MAT,
             cwd = CMD) %>%
      dplyr::select(-ID2, -Latitude, -Longitude, -Elevation)
    
    # Calculate tmax and tmin
    lgm$tmax <- rowMeans(dplyr::select(lgm, Tmax01,Tmax02,Tmax03,Tmax04,Tmax05,Tmax06,Tmax07,Tmax08,Tmax09,Tmax10,Tmax11,Tmax12))
    lgm$tmin <- rowMeans(dplyr::select(lgm, Tmin01,Tmin02,Tmin03,Tmin04,Tmin05,Tmin06,Tmin07,Tmin08,Tmin09,Tmin10,Tmin11,Tmin12))
    
    
    colnames(lgm)[-1] <- paste0(colnames(lgm)[-1], "_lgm")
    
    # Calculate bioclim variables
   bioclim_lgm = climates::bioclim2(tmin = as.data.frame(dplyr::select(lgm, Tmin01_lgm,
                                                               Tmin02_lgm,
                                                          Tmin03_lgm,
                                                          Tmin04_lgm,
                                                          Tmin05_lgm,
                                                          Tmin06_lgm,
                                                          Tmin07_lgm,
                                                          Tmin08_lgm,
                                                          Tmin09_lgm,
                                                          Tmin10_lgm,
                                                          Tmin11_lgm,
                                                          Tmin12_lgm)),
                       ## Remove site column
                       tmax = as.data.frame(dplyr::select(lgm, Tmax01_lgm,
                                                          Tmax02_lgm,
                                                          Tmax03_lgm,
                                                          Tmax04_lgm,
                                                          Tmax05_lgm,
                                                          Tmax06_lgm,
                                                          Tmax07_lgm,
                                                          Tmax08_lgm,
                                                          Tmax09_lgm,
                                                          Tmax10_lgm,
                                                          Tmax11_lgm,
                                                          Tmax12_lgm)),
                       prec = as.data.frame(dplyr::select(lgm, PPT01_lgm,
                                                          PPT02_lgm,
                                                          PPT03_lgm,
                                                          PPT04_lgm,
                                                          PPT05_lgm,
                                                          PPT06_lgm,
                                                          PPT07_lgm,
                                                          PPT08_lgm,
                                                          PPT09_lgm,
                                                          PPT10_lgm,
                                                          PPT11_lgm,
                                                          PPT12_lgm)),
                       files.as.inputs = FALSE)
   
   
   colnames(bioclim_lgm) <- paste0(colnames(bioclim_lgm), "_", "lgm")
   lgm <- cbind(lgm, bioclim_lgm)
   
   
   # Future climate - RCP85 with climateWNA
   rcp85 <- read_csv("./data/gis/climate_data/climateWNA/GBS tree climate data Climate WNA 2018_06_12_15GCM-Ensemble_rcp85_2085MSY.csv") %>%
     rename(id = ID1 , 
            tmax_sum = Tmax_sm,
            tmin_winter = Tmin_wt,
            tave = MAT,
            cwd = CMD) %>%
     dplyr::select(-ID2, -Latitude, -Longitude, -Elevation)
   
   # Calculate tmax and tmin
   rcp85$tmax <- rowMeans(dplyr::select(rcp85, Tmax01,Tmax02,Tmax03,Tmax04,Tmax05,Tmax06,Tmax07,Tmax08,Tmax09,Tmax10,Tmax11,Tmax12))
   rcp85$tmin <- rowMeans(dplyr::select(rcp85, Tmin01,Tmin02,Tmin03,Tmin04,Tmin05,Tmin06,Tmin07,Tmin08,Tmin09,Tmin10,Tmin11,Tmin12))
   
   
   colnames(rcp85)[-1] <- paste0(colnames(rcp85)[-1], "_rcp85")
   
   # Calculate bioclim variables
   bioclim_rcp85 = climates::bioclim2(tmin = as.data.frame(dplyr::select(rcp85, Tmin01_rcp85,
                                                                       Tmin02_rcp85,
                                                                       Tmin03_rcp85,
                                                                       Tmin04_rcp85,
                                                                       Tmin05_rcp85,
                                                                       Tmin06_rcp85,
                                                                       Tmin07_rcp85,
                                                                       Tmin08_rcp85,
                                                                       Tmin09_rcp85,
                                                                       Tmin10_rcp85,
                                                                       Tmin11_rcp85,
                                                                       Tmin12_rcp85)),
                                    ## Remove site column
                                    tmax = as.data.frame(dplyr::select(rcp85, Tmax01_rcp85,
                                                                       Tmax02_rcp85,
                                                                       Tmax03_rcp85,
                                                                       Tmax04_rcp85,
                                                                       Tmax05_rcp85,
                                                                       Tmax06_rcp85,
                                                                       Tmax07_rcp85,
                                                                       Tmax08_rcp85,
                                                                       Tmax09_rcp85,
                                                                       Tmax10_rcp85,
                                                                       Tmax11_rcp85,
                                                                       Tmax12_rcp85)),
                                    prec = as.data.frame(dplyr::select(rcp85, PPT01_rcp85,
                                                                       PPT02_rcp85,
                                                                       PPT03_rcp85,
                                                                       PPT04_rcp85,
                                                                       PPT05_rcp85,
                                                                       PPT06_rcp85,
                                                                       PPT07_rcp85,
                                                                       PPT08_rcp85,
                                                                       PPT09_rcp85,
                                                                       PPT10_rcp85,
                                                                       PPT11_rcp85,
                                                                       PPT12_rcp85)),
                                    files.as.inputs = FALSE)
   
   
   colnames(bioclim_rcp85) <- paste0(colnames(bioclim_rcp85), "_", "rcp85")
   rcp85 <- cbind(rcp85, bioclim_rcp85)
  
    # Mid-holocene  
    holo <- read_csv("./data/gis/climate_data/climateWNA/GBS tree climate data Climate WNA 2018_06_12_4GCM-Ensemble_midHoloceneMSY.csv") %>%
      rename(id = ID1 , 
             tmax_sum = Tmax_sm,
             tmin_winter = Tmin_wt,
             tave = MAT,
             cwd = CMD) %>%
      dplyr::select(-ID2, -Latitude, -Longitude, -Elevation)
    colnames(holo)[-1] <- paste0(colnames(holo)[-1], "_holo")
    
    
    # Join together in main climate dataset
    climate_gbs_mom <- left_join(left_join(left_join(left_join(left_join(climate_gbs_mom, 
                                                                         last1000), 
                                                                         lgm),
                                                                         holo), 
                                                                         wna_current),
                                                                         rcp85)
    
    
    dim(climate_gbs_mom)
    
    # Remove datasets no longer needed
    rm(last1000)  
    rm(lgm)
    rm(holo)
    rm(wna_current)
    rm(rcp85)
  
      
  ## Add random climate variable for testing
    climate_gbs_mom$random <- rnorm(n = nrow(climate_gbs_mom))
    climate_garden_mom$random <- rnorm(n = nrow(climate_garden_mom))
    
    
    ## Calculate average temperature
    climate_gbs_mom$tave <- (climate_gbs_mom$tmax + climate_gbs_mom$tmin) / 2 
    climate_garden_mom$tave <- (climate_garden_mom$tmax + climate_garden_mom$tmin) / 2 
    
    
    
  ## Calculating MEMs / PCNM
    
  # Calculate distance matrix
    library(geosphere)
    
    mom_lat_longs <- as.data.frame(dplyr::select(climate_gbs_mom, longitude, latitude))
    dist_mat <- matrix(NA, ncol = nrow(climate_gbs_mom), nrow = nrow(climate_gbs_mom))
    
    for(i in 1:nrow(dist_mat)){
      
      if(i %% 50 == 0){
        cat("Working on row:", i, " ... \n")
      }
      
      for(j in 1:nrow(dist_mat)){
        
        dist_mat[i, j] <- distHaversine(mom_lat_longs[i, ], 
                                        mom_lat_longs[j, ])/ 1000 # To convert to km
        
      }
    }
    
    dist_mat[1:10, 1:10]
    
  mem = vegan::pcnm(dist_mat)  
  
  # Plot out MEMs
  vegan::ordisurf(mom_lat_longs, mem$vectors[,1], bubble = 4, main = "PCNM 1")
  vegan::ordisurf(mom_lat_longs, mem$vectors[,2], bubble = 4, main = "PCNM 2")
    
  mem$values # Eigenvalues
  mem$vectors # Eigenvectors

  ## Bind positive eigenvectors to climate dataset
  climate_gbs_mom <- dplyr::bind_cols(climate_gbs_mom, as.data.frame(mem$vectors))
  
  
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
                   #   "tave_lgm",
                    #  "tmax_sum_last1000",
                   #   "tmax_sum_lgm",
                    #  "tmax_sum_holo",
                      "tmin_winter",
                    #  "tmin_winter_last1000",
                  #    "tmin_winter_lgm",
                     # "tmin_winter_holo",
                 #     "DD5",
                  #  "DD5_lgm",
                  #  "bioclim_04",
                      "random",
                      #  "tmax",
                      #"tmin",
                      # "latitude",
                      #  "longitude",
                      #  "elevation")
                      "cwd", 
                      "bioclim_04",  # Temperature seasonality
                      "bioclim_15", # Precipitation seasonality
                      "bioclim_18", # Precipitation of Warmest Quarter
                      "bioclim_19") # Precipitation of Coldest Quarter

    
    
# PCA on climate variables ------------------------------------------------

    # PCA based on climate data of maternal trees in provenance trial
      clim_pca_home <- prcomp(climate_garden_mom[,
                                                 climate_vars[-which(climate_vars == "random")]], scale = TRUE)
      summary(clim_pca_home)


    ## Format as dataframe
      clim_pca_home_vals <- as.data.frame(clim_pca_home$x)

    ## Scree plot of PCA
      fviz_screeplot(clim_pca_home)

    ## Contributions of each axis
      # fviz_contrib(clim_pca_home, choice = "var", axes = c(1))
      # fviz_contrib(clim_pca_home, choice = "var", axes = c(2))
      # fviz_contrib(clim_pca_home, choice = "var", axes = c(3))
      # fviz_contrib(clim_pca_home, choice = "var", axes = c(4))
      # fviz_contrib(clim_pca_home, choice = "var", axes = c(5))

    ## Arrows plot
      fviz_pca_var(clim_pca_home, axes = c(1,2))
      
      ## Visualize PCA results
      # fviz_contrib(clim_pca_home, choice = "var", axes = 1)
      # fviz_contrib(clim_pca_home, choice = "var", axes = 2)
      # fviz_contrib(clim_pca_home, choice = "var", axes = 3)
      # fviz_contrib(clim_pca_home, choice = "var", axes = 4)
      

      ## Save plot of contributions
      
      # fviz_pca_var(clim_pca_home, axes = c(1, 2), col.var="contrib")
      # ggsave(filename = paste0("./figs_tables/Figure S2 - climate PCA 1 ", 
      #                          Sys.Date(), ".pdf"),
      #        units = "cm",
      #        width = 20, height = 15)


    ## Predict PC values for
    # 1) Seedlings in common garden
    # 2) Chico  & IFG based on pca of climate of origin
    # 3) GBS moms
    # Need to be sure that the garden climate variables are scaled first

    dat_all_clim_scaled <- dat_all_clim
    garden_climate_scaled <- garden_climate # Scaled based on seedling climate of origin!!!
    climate_gbs_mom_scaled <- climate_gbs_mom

    for(var in climate_vars){

      dat_all_clim_scaled[, var] <- (pull(dat_all_clim_scaled, var) - mean(pull(dat_all_clim, var))) / sd(pull(dat_all_clim, var))

      garden_climate_scaled[, var] <- (pull(garden_climate_scaled, var) - mean(pull(dat_all_clim, var))) / sd(pull(dat_all_clim, var))

     climate_gbs_mom_scaled[, var] <- (pull( climate_gbs_mom_scaled, var) - mean(pull(dat_all_clim, var))) / sd(pull(dat_all_clim, var))

    }

    dat_all_clim_scaled[, climate_vars]
    garden_climate_scaled[, climate_vars]
    climate_gbs_mom_scaled[, climate_vars]

  # Predict PCA values

    dat_all_clim_pca_vals <-   as.data.frame(predict(clim_pca_home, dat_all_clim_scaled))
    dat_all_clim_pca_vals
    colnames(dat_all_clim_pca_vals) <- paste(colnames(dat_all_clim_pca_vals), "_clim", sep = "")

    garden_climate_pca_vals <-   as.data.frame(predict(clim_pca_home, garden_climate_scaled))
    garden_climate_pca_vals
    colnames(garden_climate_pca_vals) <- paste(colnames(garden_climate_pca_vals), "_clim", sep = "")

    climate_gbs_mom_pca_vals <- as.data.frame(predict(clim_pca_home, climate_gbs_mom_scaled))
    climate_gbs_mom_pca_vals
    colnames(climate_gbs_mom_pca_vals) <- paste(colnames(climate_gbs_mom_pca_vals), "_clim", sep = "")

  # Join back to main dataframes
    dat_all_clim <- bind_cols(dat_all_clim, dat_all_clim_pca_vals)
    garden_climate <- bind_cols(garden_climate, garden_climate_pca_vals)
    climate_gbs_mom <- bind_cols(climate_gbs_mom, climate_gbs_mom_pca_vals)



  ## Add PCs to list of climate vars
     climate_vars <- c(climate_vars, colnames(garden_climate_pca_vals))


# Filter down to seedlings with and without GBS data ----------------------
      
      ## Filtering down individuals who's moms don't have climate data and has GBS data
      dat_gbs_only_clim <- dplyr::filter(dat_all_clim,
                                         accession %in% climate_gbs_mom$accession)
      dim(dat_gbs_only_clim)
      
      
      # Select a 1/2 subset of adults for taining and testing
      # set.seed(1)
      # training_moms <- sample(na.omit(climate_gbs_mom$accession), 
      #                         size = length(na.omit(climate_gbs_mom$accession))/2)
      # testing_moms <- na.omit(climate_gbs_mom$accession)[!na.omit(climate_gbs_mom$accession) %in% training_moms]
      # 
      # any(training_moms %in% testing_moms) # Should both be false
      # any(testing_moms %in% training_moms)
      # 
      # dat_gbs_only_clim <- dat_gbs_only_clim %>%
      #   dplyr::filter(accession %in% training_moms)
      
    
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
  


 # pairs.panels(climate_garden_mom[, c("latitude", "longitude", 
 #                                     #"mem1", "mem2", "mem3", 
 #                                     climate_vars)])
 # 
 # pairs.panels(climate_garden_mom[, c("tmax_sum", "tmin_winter", "DD5", "tmax_sum_lgm","tmin_winter_lgm", "DD5_lgm")], ellipses = F)
 # 
 # pairs.panels(climate_garden_mom[, c("tmax_sum",
 #                                     "tmax_sum_lgm",
 #                                     "tmin_winter", 
 #                                     "tmin_winter_lgm",
 #                                     "DD5", 
 #                                     "DD5_lgm")], ellipses = F, cex = 2)
 # 
 # pairs.panels(climate_garden_mom[, c("tmax_sum",
 #                                     "tmax",
 #                                     "DD5",
 #                                     "tave",
 #                                     "bioclim_04")])
 
 #pairs.panels(dat_all_scaled[, climate_vars_dif])
 



