

# Load libraries ----------------------------------------------------------

source("./01_clean-process-explore-data.R")
source("./03_adding-climate-data.R")

dim(dat_all_scaled)

# install_github("jdstorey/qvalue")
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
  
  future_df <-  as.data.frame(cbind(ccsm4_85,
                                    cnrm_85, 
                                    ipsl_85, 
                                    miroc_85,
                                    fgoals_85,
                                    miroc_45,
                                    mpi_45))
  
  pairs.panels(future_df)
  
  ## Density plots of different scenarios for provenances
  ## Miroc much higher than the others
  future_df %>%
    gather(key = "projection", "tmax_sum") %>%
    ggplot(., aes(tmax_sum, fill = projection)) + geom_density(alpha = 0.5) + 
    theme_bw(15) 
  
  # Calculate mean increases
  future_85_mean <- mean(apply(dplyr::select(future_df, contains("85")), MARGIN = 1, mean))
  future_45_mean <- mean(apply(dplyr::select(future_df, contains("45")), MARGIN = 1, mean))
  


    

# Run FULL models with all seedlings (no genomic data) --------------------

# Convert factors
  dat_all_scaled$section <- factor(dat_all_scaled$section)
  dat_all_scaled$section_block <- factor(dat_all_scaled$section_block) 
  dat_all_scaled$accession <- factor(dat_all_scaled$accession)
  dat_all_scaled$site <- factor(dat_all_scaled$site)
  
  var = "tmax_sum_dif"


# Set formula for gam
  form <- formula(paste0("rgr ~ section_block + s(height_2014, bs =\"cr\") + s(", var, " , bs = \"cr\") + s(accession, bs = \"re\")"))
  
  
  # With all 5,000+ seedlings
  gam_all <- bam(formula = form,
                 data = dat_all_scaled,
              #   discrete = TRUE, 
                 nthreads = 8,
                 method = "fREML", family = "gaussian", 
                 control = list(trace = FALSE))
  
  summary(gam_all)

  # Plot overall model fit
  test_fit <- dat_all_scaled
  test_fit$pred <- gam_all$fitted.values
  

  ggplot(test_fit, aes(x = pred, y = rgr, bg = section_block)) + geom_point(alpha = 0.75, pch = 21) + theme_bw(15) +
    geom_abline(slope = 1, intercept = 0, lwd = 1.5, col = "forestgreen")
  
  visreg(gam_all, partial = TRUE, ylab = "5 yr height (cm)")
  
  visreg(gam_all, partial = FALSE, ylab = "5 yr height (cm)")
  
  gam.check(gam_all)
  
  ## Save gam summary to file
    # sink(file = paste0("./figs_tables/Table 1 - full_gam_model_summary_",
    #                    Sys.Date(), ".txt"))
    # summary(gam_all)
    # anova(gam_all)
    # sink()
  

#### Prediction plot of changes in height with tmax transfer  
  
  
  # Make predictions
  v <- visreg(gam_all, xvar = "tmax_sum_dif", 
              scale = "response", rug = FALSE, ylab = "Relative growth rate", plot = FALSE,
              # Set prediction values
              cond = list(section_block = "IFG_1", height_2014 = 0, accession = 1)) 
  
    # Transform x axis + assign variables for lo and hi CI
    v$fit <- v$fit %>%
      # Back transform X axis
      mutate(tmax_sum_dif = back_transform(x = tmax_sum_dif,
                                      var = "tmax_sum_dif",
                                      means = scaled_var_means_all,
                                      sds = scaled_var_sds_all))
  ## Main GGPLOT
  gg <- 
    ggplot(v$fit, aes(x = tmax_sum_dif, y = rgr)) +
    geom_ribbon(data = v$fit, aes(ymin = visregLwr, ymax = visregUpr), 
                fill = "grey80", alpha = 0.75) +
    geom_vline(aes(xintercept = 0), lty = 2) +
  #  geom_vline(xintercept = future_45_mean) + 
    geom_vline(xintercept = future_85_mean, lwd = 1.5) +
    geom_line(data = v$fit, aes(x = tmax_sum_dif, y = visregFit), lwd = 2,
              col = "forestgreen") +
    scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5)) +
   # ylab(expression(Relative~growth~rate~(cm~cm^-1~yr^-1))) +
    ylab("Relative growth rate") +
    xlab("Tmax transfer distance") +
   # ylim(c(110, 200)) +
    theme_bw(15) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.margin = (margin(1.5,1.5,1.5,1.5, "cm"))) +
    NULL
  
  gg
  
  
# Predicting changes in height based on degree increase -------------------

  # Set degree increases to compare for no change - medium & high emissions
  degrees <- c(0, future_45_mean, future_85_mean)
  
  newdata <- data.frame(section_block = v$fit$section_block[1],
                        height_2014 = v$fit$height_2014[1],
                        accession = v$fit$accession[1],
                        tmax_sum_dif_unscaled = degrees,
                        tmax_sum_dif = forward_transform(x = degrees,
                                                         var = "tmax_sum_dif",
                                                         means = scaled_var_means_all,
                                                         sds = scaled_var_sds_all))
  
  newdata
  
  newdata$pred <- predict(gam_all, newdata = newdata, se.fit = TRUE)$fit
  newdata$se <- predict(gam_all, newdata = newdata, se.fit = TRUE)$se.fit
  
  newdata$lwr <- newdata$pred - newdata$se * 1.96
  newdata$upr <- newdata$pred + newdata$se * 1.96
  
  newdata
  
## % change medium emissions scenario
  
  ## Mean
  (newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[2]] - 
    newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[1]]) / newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[1]] * 100
  
  ## Lwr
  (newdata$lwr[newdata$tmax_sum_dif_unscaled == degrees[2]] - 
      newdata$lwr[newdata$tmax_sum_dif_unscaled == degrees[1]]) / newdata$lwr[newdata$tmax_sum_dif_unscaled == degrees[1]] * 100
  
  ## Upr
  (newdata$upr[newdata$tmax_sum_dif_unscaled == degrees[2]] - 
      newdata$upr[newdata$tmax_sum_dif_unscaled == degrees[1]]) / newdata$upr[newdata$tmax_sum_dif_unscaled == degrees[1]] * 100
  
  
  ## % change in 5 degree increase
  
  ## Mean
  (newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[3]] - 
      newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[1]]) / newdata$pred[newdata$tmax_sum_dif_unscaled == degrees[1]] * 100
  
  
  
  

# Make spatial predictions across landscape -------------------------------

  # Load libraries ----------------------------------------------------------
  
  library(raster)
  library(rasterVis)
  
  
# Loading in rasters ------------------------------------------------------
  
  ## Read in elevation dem and create a hillshade for mapping
  dem <- raster("./data/gis/dem/ClimateNA_DEM_cropped.tif")
  
  slope = raster::terrain(dem, opt='slope')
  aspect = raster::terrain(dem, opt='aspect')
  hill = hillShade(slope, aspect, 40, 270)
  
  ## Load in california outline
  cali_outline <- readShapePoly("./data/gis/california_outline/california_outline.shp",
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  lobata_range <- readShapePoly("./data/gis/valley_oak_range/qlobata_refined_grouped.shp",
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  lobata_range_rough <- readShapePoly("./data/gis/valley_oak_range/querloba.shp",
                                      proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  ## Creates vector with list of file paths to all .tif raster files
   ## Directory path where future climate scenarios are located
  dir_name <- "./data/gis/climate_data/BCM/future/"
  raster_files <- list.files(dir_name, full.names = TRUE, recursive = TRUE)
  raster_files <- raster_files[grep("*[[:digit:]].tif$", raster_files)] # Only those with .tif extension
  raster_files <- raster_files[grep("jja_ave", raster_files)] # Only tmax sum files
  
  raster_files
  
  future_stack <- stack()
  
  ## Loop through raster files and add climate values to dat_mom dataframe
  for(file in raster_files){
    
    cat("Working on file:", file, "...\n")
    ## Load in raster
    future_stack <- raster::stack(future_stack, raster(file))
  
  }
  
  future_stack
  
  ## Clean up names
  names(future_stack) <-  unlist(lapply(strsplit(x = gsub(pattern = "tmx2070_2099jja_ave_", 
                                                          replacement = "", names(future_stack)),
                                                 split = "_1"), "[", 1))
  
 
  
  tmax_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmx1951_1980jja_ave_HST_1513103038/tmx1951_1980jja_ave_HST_1513103038.tif")
  
  
# Calculate difference between future and current
  future_stack_tmax_dif <- future_stack - tmax_rast
  names(future_stack_tmax_dif) <- names(future_stack) # Reassign names
  
  # Scale values
  future_stack_tmax_dif_scaled <- (future_stack_tmax_dif - scaled_var_means_all["tmax_sum_dif"]) /
    scaled_var_sds_all["tmax_sum_dif"]
  
  
  # Make data frame for prediction
  newdata_rast <- data.frame(section_block = v$fit$section_block[1],
                        height_2014 = v$fit$height_2014[1],
                        accession = v$fit$accession[1])
  
  ## Get null prediction of height with no change in tmax
    future_null <- tmax_rast
    names(future_null) <- "tmax_sum_dif"
    
    # Need to scale tmax when actual transfer distance is 0
    future_null[!is.na(future_null)] <- forward_transform(0, var = "tmax_sum_dif",
                                                          means = scaled_var_means_all,
                                                          sds = scaled_var_sds_all)
    future_null
    
    future_null_height <- predict(future_null,
                                  model = gam_all,
                                  const = newdata_rast,
                                  progress = "text")
    
    future_null_height
    summary(future_null_height)
  
    
  # Initialize stack for height changes
    future_stack_height_change <- future_stack_tmax_dif
  
  # Loop over scenarios
  for(scenario in 1:nlayers(future_stack_tmax_dif)){
    
    cat("Working on scenario:", names(future_stack_tmax_dif)[scenario], "... \n")
    
    ## Use scaled raster
    rast_temp <- future_stack_tmax_dif_scaled[[scenario]]
    names(rast_temp) <- "tmax_sum_dif" # Rename so predict function works
    
    rast_temp_height <- predict(rast_temp,
                        model = gam_all,
                        const = newdata_rast,
                        progress = "text")
    
    ## Calculate % change in height
    rast_temp_height <- (rast_temp_height - future_null_height) / future_null_height * 100
    
    future_stack_height_change[[scenario]] <- rast_temp_height
    names(future_stack_height_change)[scenario] <- names(future_stack_tmax_dif)[scenario] # Reassign name
    
  } # End scenario loop
    
    
  # Plot summaries across scenarios  
   # histogram(future_stack_height_change)
  # bwplot(future_stack_height_change)
  
  ## Write Raster to file to avoid having to run it again
    # writeRaster(x = future_stack_height_change,
    #             filename = paste0("./output/future_stack_height_change_", Sys.Date(), ".tif"),
    #             progress = "text")
    
    
  ## Average across low and high emissions scenarios
    future_stack_height_change_85 <- mean(future_stack_height_change[[c("CCSM4_rcp85", 
                                                                        "CNRM_rcp85",
                                                                        "Fgoals_rcp85",
                                                                        "IPSL_rcp85",
                                                                        "MIROC_rcp85")]])
    
    summary(future_stack_height_change_85 )
    
    future_stack_height_change_45 <- mean(future_stack_height_change[[c("MPI_rcp45", 
                                                                        "MIROC_rcp45")]])
    
    summary(future_stack_height_change_45)
    
    ## Project to latlong
    future_stack_height_change_85 <- projectRaster(future_stack_height_change_85, 
                                                   crs = CRS("+proj=longlat +datum=WGS84"))
    names(future_stack_height_change_85) <- "RCP85"
    
    future_stack_height_change_45 <- projectRaster(future_stack_height_change_45, 
                                                   crs = CRS("+proj=longlat +datum=WGS84"))
    names(future_stack_height_change_45) <- "RCP45"
    
    ## Stack them
    future_stack_height_change_45_85 <-  stack(#future_stack_height_change_45,
          future_stack_height_change_85)
    
    ## Crop them
    future_stack_height_change_45_85 <- mask(future_stack_height_change_45_85, cali_outline)
    
    
    
  ## Options for levelplot
    pixel_num = 1e6 ## Can make resolution better by making this 1e6 
    
    myTheme <- rasterTheme(region = brewer.pal('RdYlBu', n = 9))
    
   breaks <- c(seq(-8, -4, by = 0.2))
    
  ## Plots of changes in height  
    levelplot(future_stack_height_change_45_85,
              margin = FALSE,
              maxpixels = pixel_num,
              par.settings = myTheme, 
              at = breaks) + 
    latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey10")) + 
    latticeExtra::layer(sp.polygons(lobata_range, col = "grey50")) + 
    latticeExtra::layer(sp.polygons(lobata_range_rough, col = "black", lwd = 2.75))
  
    
    # Get values by region
    height_change_by_region <- raster::extract(future_stack_height_change_45_85,
                                              y = lobata_range, df = TRUE,
                                              fun = mean, na.rm = TRUE)
  
    height_change_by_region$region <- lobata_range$REGION
    
    height_change_by_region <- height_change_by_region %>%
                            gather(key = "scenario", 
                                   value = "height_change", 
                                   -ID, -region)
    
  # Barplot
 
    height_change_by_region <- raster::extract(future_stack_height_change_45_85,
                                               y = lobata_range, df = TRUE)
    height_change_by_region$ID <- factor(height_change_by_region$ID)
    
    ggplot(height_change_by_region, aes(ID, RCP85, group = ID, bg = ID)) + 
      geom_boxplot() + theme_bw(30) + 
      scale_x_discrete(labels = lobata_range$REGION) + 
      scale_fill_discrete(guide=FALSE) +
      xlab("Region") + ylab("Change in RGR")

    
  
  
# Graveyard ---------------------------------------------------------------
  
  
  # Raster processing -------------------------------------------------------
  
  
  ## Read in elevation dem and create a hillshade for mapping
  dem <- raster("./data/gis/dem/ClimateNA_DEM_cropped.tif")
  
  slope = raster::terrain(dem, opt='slope')
  aspect = raster::terrain(dem, opt='aspect')
  hill = hillShade(slope, aspect, 40, 270)
  
  
  ## Get general file names without extension
  general_files <- c(list.files(path="./data/gis/climate_data/NA_NORM_8110_Bioclim_netCDF/",
                                pattern = "_cropped.tif$", full.names = FALSE),
                     list.files(path="./data/gis/climate_data/NA_ENSEMBLE_rcp85_2080s_Monthly_netCDF/",
                                pattern = "_cropped.tif$", full.names = FALSE))
  general_files <- gsub(pattern = ".tif", replacement = "", general_files)
  
  var_names <- paste(gsub(pattern = "_cropped", "", general_files), "_dif", sep = "")
  
  general_files <- paste("/", general_files, sep = "") ## Need this for grep to be specific
  
  ## Read full files names of current climates
  current_files <- c(list.files(path="./data/gis/climate_data/NA_NORM_8110_Bioclim_netCDF/",
                                pattern = "_cropped.tif$", full.names = TRUE),
                     list.files(path="./data/gis/climate_data/NA_NORM_8110_Monthly_netCDF/",
                                pattern = "_cropped.tif$", full.names = TRUE))
  
  ## Read full file names of future climate
  future_files <- c(list.files(path="./data/gis/climate_data/NA_ENSEMBLE_rcp85_2080s_Bioclim_netCDF/",
                               pattern = "_cropped.tif$", full.names = TRUE),
                    list.files(path="./data/gis/climate_data/NA_ENSEMBLE_rcp85_2080s_Monthly_netCDF/",
                               pattern = "_cropped.tif$", full.names = TRUE))
  
  
  
  # Calculate Tmax and Tmin -------------------------------------------------
  
  
  ### Calculate max and min temperatature of the warmest month
  dat_all <- dat_all %>%
    mutate(Tmax = apply( dplyr::select(., Tmax01:Tmax12), 1, max))
  
  
  ### Tmax
  Tmax_current_rast <- stack()
  Tmax_future_rast <- stack()
  
  for(var in general_files[grep("Tmax", general_files)]){
    
    current_path <- current_files[grep(var, current_files)]
    cat("Reading in as current file:", current_path, "... \n")
    
    ## Load in current raster
    current_rast <- raster(current_path) 
    
    future_path <- future_files[grep(var, future_files)]
    cat("Reading in as future file:", future_path, "... \n")
    
    ## Load in future raster 
    future_rast <-  raster(future_path) 
    
    ### Save to respective stack
    Tmax_current_rast <- raster::stack(Tmax_current_rast, current_rast)
    Tmax_future_rast <- raster::stack(Tmax_future_rast, future_rast)
    
  }
  
  ## Calculate difference between future and current values
  
  Tmax_current_rast_max <- max(Tmax_current_rast)
  names(Tmax_current_rast_max) <- "Tmax"
  
  Tmax_dif_rast <- max(Tmax_future_rast) - Tmax_current_rast_max
  names(Tmax_dif_rast) <- "Tmax_dif"
  plot(Tmax_dif_rast)
  
  
  ### Tmin
  Tmin_current_rast <- stack()
  Tmin_future_rast <- stack()
  
  for(var in general_files[grep("Tmin", general_files)]){
    
    current_path <- current_files[grep(var, current_files)]
    cat("Reading in as current file:", current_path, "... \n")
    
    ## Load in current raster
    current_rast <- raster(current_path) 
    
    future_path <- future_files[grep(var, future_files)]
    cat("Reading in as future file:", future_path, "... \n")
    
    ## Load in future raster 
    future_rast <-  raster(future_path) 
    
    ### Save to respective stack
    Tmin_current_rast <- raster::stack(Tmin_current_rast, current_rast)
    Tmin_future_rast <- raster::stack(Tmin_future_rast, future_rast)
    
  }
  
  ## Calculate difference between future and current values
  
  Tmin_current_rast_min <- min(Tmin_current_rast)
  names(Tmin_current_rast_min) <- "Tmin"
  
  Tmin_dif_rast <- min(Tmin_future_rast) - Tmin_current_rast_min
  names(Tmin_dif_rast) <- "Tmin_dif"
  
  plot(Tmin_dif_rast) 
  
  
  
  
  ## Loop through climate variables and calculate difference between future 
  ## and current climate - Tmax and Tmin addded in later
  
  
  
  ## Initialize stack of rasters  
  dif_rast <- stack()
  cur_rast <- stack()
  x = 1 # For looping through var names in plot title
  
  ## Only loop through climate variables that are used in the models to speed up process
  for(var in general_files[apply(sapply(paste(climate_vars, "_", sep =""), 
                                        function(x) grepl(x, general_files)), 1, any)]){
    
    cat("\n ******* \n")
    cat("Calculating difference in:", var, "... \n")
    
    current_path <- current_files[grep(var, current_files)]
    cat("Reading in as current file:", current_path, "... \n")
    
    ## Load in current raster
    current_rast <- raster(current_path) 
    
    future_path <- future_files[grep(var, future_files)]
    cat("Reading in as future file:", future_path, "... \n")
    
    ## Load in future raster 
    future_rast <-  raster(future_path) 
    
    ### Calculate difference as future minus current
    
    difference_rast <- future_rast - current_rast
    
    names(difference_rast) <- var
    
    ## Plot output
    # plot(difference_rast, main = paste("Difference in:", var_names[x]))
    
    ## Add to stack of rasters
    dif_rast <- raster::stack(dif_rast, difference_rast)
    cur_rast <- raster::stack(cur_rast, current_rast)
    
    x = x+1
  }  
  
  names(dif_rast) <- gsub("cropped", "dif", gsub("^X.", "", names(dif_rast)))
  names(cur_rast) <- gsub("_dif", "", names(dif_rast))
  
  ### Add in Tmax and Tmin difference rasts
  dif_rast <- raster::stack(dif_rast, Tmax_dif_rast)
  dif_rast <- raster::stack(dif_rast, Tmin_dif_rast)
  
  ## Add Tmax and Tmin to current rasters
  cur_rast <- raster::stack(cur_rast, Tmax_current_rast_max)
  cur_rast <- raster::stack(cur_rast, Tmin_current_rast_min)
  
  
  dif_rast_raw  <- dif_rast # Intermediate object so don't have to reload if change masking
  
  
  
  # Cropping based on lobata range ------------------------------------------
  
  dif_rast = dif_rast_raw
  # dif_rast <- raster::mask(dif_rast_raw, rgeos::gBuffer(lobata_range,                                                      width = 0, byid = TRUE))
  ## If including current conditions
  # dif_rast <- raster::stack(dif_rast_raw, cur_rast)
  
  
  
  # Scaling elevation data --------------------------------------------------
  
  elev_rast <- dem
  names(elev_rast) <- "elev"
  elev_rast <- log(elev_rast + 10) ## Similar to how it was handled in processing script
  
  elev_df <- data.frame(elev = raster::extract(elev_rast, 1:ncell(elev_rast)))
  summary(elev_df)
  
  
  # Extract and scale values ------------------------------------------------
  
  ## Extract values from raster into dataframe
  extracted_vals = raster::extract(dif_rast, 1:ncell(dif_rast))
  extracted_vals <- as.data.frame(extracted_vals)
  
  ## Add in elevation
  extracted_vals <- cbind(extracted_vals, elev_df)
  
  ### Scale extracted values to match with scaled variables
  extracted_vals_scaled <- extracted_vals
  
  for(var in c(climate_vars_dif, "elev")){
    extracted_vals_scaled[, var] <- (extracted_vals[, var] - scaled_var_means[var]) / 
      scaled_var_sds[var]
  }
  
  summary(extracted_vals_scaled)
  
  
  
  # Histogram of predicted climate dif vs garden dif ------------------------
  
  # # Make histogram of predicted climate differences vs. those from the gardens
  #   extracted_vals_nona <- extracted_vals[!is.na(extracted_vals[, 1]), ] # Remove NAs
  # 
  #   dif_df <- tidyr::gather(extracted_vals_nona, key = "climate_var_dif",
  #                         value = "dif")
  #   dif_df$type <-  "Predicted climate"
  # 
  #   observed <- tidyr::gather(dplyr::select(dat_all, climate_vars_dif),
  #                             key = "climate_var_dif",
  #                         value = "dif")
  #   observed$type = "Provenance differences"
  # 
  #   dif_df <- rbind(dif_df, observed)
  # 
  #   for(var in levels(factor(dif_df$climate_var_dif))){
  # 
  #     p = ggplot(dif_df[dif_df$climate_var_dif == var, ], aes(x = dif, fill = type)) +
  #         geom_density(col = "black", alpha = 0.5) + ggtitle(var) + xlab(var) +
  #         theme_bw() + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
  # 
  #     print(p)
  # 
  #   }
  # 
  # 
  # ## Plots of predicted changes
  # 
  # ## CMD
  # levelplot(dif_rast[["CMD_dif"]], contour = FALSE, margin = FALSE, main = "CMD_dif",
  #           par.settings = viridisTheme()) +  layer(sp.polygons(cali_outline, lwd=2))  +
  #   layer(sp.polygons(lobata_range, fill = "black", alpha = 0.35))
  # 
  # ## TMin
  # levelplot(dif_rast[["Tmin_dif"]], contour = FALSE, margin = FALSE, main = "Tmin_dif",
  #           par.settings = viridisTheme()) 
  # 
  # ## Tmax
  # rasterVis::levelplot(dif_rast[["Tmax_dif"]], contour = FALSE, margin = FALSE, main = "Tmax_dif", par.settings = viridisTheme())
  # 
  # ## DD5
  # levelplot(dif_rast[["DD5_dif"]], contour = FALSE, margin = FALSE, main = "DD5_dif",
  #           par.settings = viridisTheme()) +  layer(sp.polygons(cali_outline, lwd=2))  +
  #   layer(sp.polygons(lobata_range, fill = "black", alpha = 0.35))
  
  
  
  

  
  
  
  
  
  
  
  
# Graveyard ---------------------------------------------------------------

  # Function for plotting predictions ---------------------------------------
  
  plot_pred <- function(xvar,
                        xlab = "",
                        add_points = TRUE,
                        add_residuals = FALSE,
                        ylim,
                        save = FALSE,
                        path = "./",
                        filename = "",
                        main = ""){
    
    # Test to see what type of variable it is - factor or numeric
    type <- ifelse(is.factor(pull(dat_all_scaled, xvar)), "factor", "numeric")
    
    # Make predictions
    v <- visreg(gam_all, xvar = xvar, 
                scale = "response", rug = FALSE, ylab = "5 yr height (cm)", plot = FALSE,
                cond = list(section_block = "IFG_1", height_2014 = 0, accession = 1))
    
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
      xlab(xlab) +
      ylim(ylim) +
      theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
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
        geom_ribbon(data = v$fit, aes(ymin = visregLwr, ymax = visregUpr), fill = "grey80", alpha = 0.75) +
        geom_line(data = v$fit, aes(x = get(xvar), y = visregFit), lwd = 2, col = "forestgreen") 
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
                     text = element_text(size = 15))
    
    # Printing to file
    if(save){
      ggsave(paste0(path, filename))
    }
    
    print(gg)
    
    return(gg)
    
  } # End plot_pred function
  
  

  