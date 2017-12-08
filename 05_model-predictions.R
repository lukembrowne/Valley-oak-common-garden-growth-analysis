

# Load libraries ----------------------------------------------------------

library(raster)
library(rasterVis)


# Raster processing -------------------------------------------------------


  ## Read in elevation dem and create a hillshade
    dem <- raster("./data/gis/dem/ClimateNA_DEM_cropped.tif")
    
    slope = terrain(dem, opt='slope')
    aspect = terrain(dem, opt='aspect')
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

  Tmax_dif_rast <- max(Tmax_future_rast) - max(Tmax_current_rast)
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
  
  Tmin_dif_rast <- min(Tmin_future_rast) - min(Tmin_current_rast)
  names(Tmin_dif_rast) <- "Tmin_dif"
  
  plot(Tmin_dif_rast) 


## Loop through climate variables and calculate difference between future and current climate
  
  

  ## Initialize stack of rasters  
  dif_rast <- stack()
  cur_rast <- stack()
  x = 1 # For looping through var names in plot title
  
  ## Only loop through climate variables that are used in the models to speed up process
  ## Tmax and Tmin added in later
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
  
  names(dif_rast) <- gsub("cropped", "dif", gsub("X.", "", names(dif_rast)))
  names(cur_rast) <- gsub("_dif", "", names(dif_rast))
  
  ### Add in Tmax and Tmin difference rasts
  dif_rast <- raster::stack(dif_rast, Tmax_dif_rast)
  dif_rast <- raster::stack(dif_rast, Tmin_dif_rast)
  
  
  dif_rast_raw  <- dif_rast # Intermediate object so don't have to reload if change masking



# Cropping based on lobata range ------------------------------------------

  dif_rast = dif_rast_raw
  # dif_rast <- raster::mask(dif_rast_raw, rgeos::gBuffer(lobata_range,                                                      width = 0, byid = TRUE))
  

# Extract and scale values ------------------------------------------------

  ## Extract values from raster into dataframe
  extracted_vals = raster::extract(dif_rast, 1:ncell(dif_rast))
  extracted_vals <- as.data.frame(extracted_vals)


  ### Scale extracted values to match with scaled variables

  extracted_vals_scaled <- extracted_vals
  
    for(var in climate_vars_dif){
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



# Make layers for PCA model -----------------------------------------------

  ### PCA model  
  ### Make new rasters of each PC by calculating their loadings
  # predictions <- predict(clim_pca, newdata = extracted_vals_scaled)
  # predictions <- as.data.frame(predictions)
  # 
  # summary(predictions)
  # 
  # pc1_rast <- dif_rast[[1]]
  # 
  # values(pc1_rast) <- predictions$PC1
  # 
  # levelplot(pc1_rast, contour = TRUE, margin = FALSE)
  # 


# Model predictions -------------------------------------------------------

  
plot_future_growth <- function(gam_mod){
  
  ## Make prediction dataframe
  predictions <- extracted_vals_scaled[, climate_vars_dif] # Doesn't matter for univariate models if there are multiple climatic variables in this DF
  summary(predictions)
  
  predictions$height_2016 <-  0 ## Use average seedling height
  predictions$cluster_assigned <- 1 ## If using genetic model
  
  ## Make a null predictions set for comparison
  predictions_null <- predictions
  
  for(var in climate_vars_dif){
    ## Set predictions at scaled value that equals 0 change in climate
    predictions_null[, var] <- ( 0 - scaled_var_means[[var]]) / scaled_var_sds[[var]]
  }
  summary(predictions_null)
  
  
  ## Generate predictions
  height_predictions =  predict(gam_mod, predictions, 
                                re.form = NA, type = "response",
                                se.fit = FALSE)
  
  summary(height_predictions)
  
  height_predictions_null <-  predict(gam_mod, predictions_null, 
                                      re.form = NA, type = "response",
                                      se.fit = FALSE)
  summary(height_predictions_null)
  
  ## Save values into raster
  height_rast <- dif_rast[[1]] ## Initialize raster
  height_rast_null <- dif_rast[[1]]
  
  values(height_rast) <- NA # Just in case the next line doesn't work
  values(height_rast) <- c(height_predictions) # Save in real predictions
  
  values(height_rast_null) <- NA # Just in case the next line doesn't work
  values(height_rast_null) <- c(height_predictions_null) # Save in real predictions
  
  ## Percent change in growth rates in comparison to no change 
  height_rast_change <- ((height_rast - height_rast_null) / height_rast ) * 100
  
  hist(height_rast_change)
  
  
  ## Plot of predicted growth rates
  hsTheme <- modifyList(GrTheme(), list(regions=list(alpha=.15)))
  heightTheme <- modifyList(RdBuTheme(), list(regions=list(alpha=1)))
  
  pixel_num = 1e5 ## Can make resolution better by making this 1e6 
  
  ## Plot
  cat("Plotting % change in growth rates...")
  
  p = levelplot(height_rast_change, 
    #  levelplot(raster::mask(height_rast, rgeos::gBuffer(lobata_range, width = 0, byid = TRUE)),
            contour = FALSE, margin = FALSE, par.settings = heightTheme, 
            maxpixels = pixel_num, colorkey = TRUE, main = "% Change in Relative Growth Rate") +
    levelplot(hill,  margin = FALSE, par.settings=hsTheme, maxpixels = pixel_num) +
    latticeExtra::layer(sp.polygons(cali_outline, lwd=2)) +
    latticeExtra::layer(sp.polygons(lobata_range,fill = "black", alpha = 0.15))
  
  print(p)
  
  
  ## Histogram of changes in relative height
  
  ## Across all of california
  
  height_change_range = unlist(raster::extract(height_rast_change, lobata_range))
  height_change_cali = unlist(raster::extract(height_rast_change, 1:ncell(height_rast)))
  
  summary(height_change_range)
  summary(height_change_cali)
  
  ## Growth within range
  p = ggplot(data.frame(height_change_range),
         aes(height_change_range)) + 
    geom_histogram(color = "black", fill = "steelblue2") + theme_bw() +
    geom_vline(xintercept = 0, lty = 2) + 
    xlab("% Change in relative growth rate within range") + 
    theme(plot.margin = margin(1, 1, 1, 1, "cm"))
  
  print(p)
  
  ## Growth within Cali
  p = ggplot(data.frame(height_change_cali),
         aes(height_change_cali)) + 
    geom_histogram(color = "black", fill = "steelblue2") + theme_bw() +
    geom_vline(xintercept = 0, lty = 2) +
    xlab("% Change in relative growth rate within California") + 
    theme(plot.margin = margin(1, 1, 1, 1, "cm"))
  
  print(p)
  
} # End function
  
  ## CMD
  plot_future_growth(fit_CMD$gam)
  
  ## Tmax
  plot_future_growth(fit_Tmax$gam)
  
  ## Tmin
  plot_future_growth(fit_Tmin$gam)
  
  ## DD5
  plot_future_growth(fit_DD5$gam)
  
  

## HEIGHT PREDICTIONS - Use PC scores to make predictions about height




