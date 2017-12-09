

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



# Model predictions -------------------------------------------------------

  
plot_future_growth <- function(gam_mod, PCA = FALSE){
  
  ## Make prediction dataframe
  # Doesn't matter for univariate models if there are >1 climatic variables in this DF
  if(PCA == FALSE){
    predictions <- extracted_vals_scaled[, c(climate_vars_dif, "elev")] 
     summary(predictions) 
  }  else if(PCA == TRUE){
    ## PCA model
    ## Make new rasters of each PC by calculating their loadings
    
    predictions <- predict(clim_pca, newdata = extracted_vals_scaled)
    predictions <- as.data.frame(predictions)
    
    colnames(predictions) <- paste(colnames(predictions), "_clim_dif", sep = "")
    
    summary(predictions)
  }
  
  predictions$height_2016 <-  0 ## Use average seedling height
  predictions$cluster_assigned <- 1 ## If using genetic model
  
  ## Make a null predictions set for comparison
  predictions_null <- predictions
  
  for(var in climate_vars_dif){
    ## Set predictions at scaled value that equals 0 change in climate
    predictions_null[, var] <- ( 0 - scaled_var_means[[var]]) / scaled_var_sds[[var]]
  }
  summary(predictions_null)
  
  ## Convert to PCs
  if(PCA == TRUE){
    predictions_null <- predict(clim_pca, newdata = predictions_null)
    
    colnames(predictions_null) <- paste(colnames(predictions_null), "_clim_dif", sep = "")
    
    predictions_null <- as.data.frame(predictions_null)
    
    predictions_null$height_2016 <-  0 ## Use average seedling height
  }
  
  ## Generate predictions
  rgr_predictions =  predict(gam_mod, predictions, 
                                re.form = NA, type = "response",
                                se.fit = TRUE)
  
  summary(rgr_predictions$fit)
  hist(rgr_predictions$fit)
  
  ## Null model
  rgr_predictions_null <-  predict(gam_mod, predictions_null, 
                                      re.form = NA, type = "response",
                                      se.fit = TRUE)
  summary(rgr_predictions_null$fit)
  hist(rgr_predictions_null$fit)
  
  ## Calculating points where SEs don't overlap
  rgr_predictions$upper <- (rgr_predictions$fit + 2 * rgr_predictions$se.fit)
  rgr_predictions_null$lower <- (rgr_predictions_null$fit - 2 *
                                   rgr_predictions_null$se.fit)
  rgr_predictions$noSE_overlap <- rgr_predictions$upper <= rgr_predictions_null$lower
  
  cat("Whether SE overlap or not - if TRUE  == SE doesn't overlap \n")
  print(table(rgr_predictions$noSE_overlap)) ## IF TRUE THIS MEANS SE DON'T OVERLAP
  

  ## Save values into raster
  rgr_rast <- dif_rast[[1]] ## Initialize raster
  rgr_rast_null <- dif_rast[[1]]
  rgr_rast_SEoverlap <- dif_rast[[1]]
  
  values(rgr_rast) <- NA # Just in case the next line doesn't work
  values(rgr_rast) <- c(rgr_predictions$fit) # Save in real predictions
  
  values(rgr_rast_null) <- NA # Just in case the next line doesn't work
  values(rgr_rast_null) <- c(rgr_predictions_null$fit) # Save in real predictions
  
  values(rgr_rast_SEoverlap) <- NA
  values(rgr_rast_SEoverlap) <- c(rgr_predictions$noSE_overlap)
  
  ## Percent change in growth rates in comparison to no change 
  ## Formula == (new value - old value) / old value
  rgr_rast_change <- ((rgr_rast - rgr_rast_null) / rgr_rast_null ) * 100
  
  ## Getting rid of outliers
  rgr_rast_change[rgr_rast_change < quantile(rgr_rast_change, 0.025) | 
                    rgr_rast_change > quantile(rgr_rast_change, 0.975)] <- NA
  
  rgr_rast_change[rgr_rast_change < -100 | 
                    rgr_rast_change > 100] <- NA
  summary(rgr_rast_change)
  hist(rgr_rast_change)
  
  ## Problem cells [1] 1068038 1361177 1486313
  
  
  ## Plot of predicted growth rates
  hsTheme <- modifyList(GrTheme(), list(regions=list(alpha=.15)))
  rgrTheme <- modifyList(RdBuTheme(), list(regions=list(alpha=1)))
  
  pixel_num = 1e6 ## Can make resolution better by making this 1e6 
  
  ## Plot
  cat("Plotting % change in growth rates...")
  
  p = levelplot(rgr_rast_change, 
    #  levelplot(raster::mask(rgr_rast, rgeos::gBuffer(lobata_range, width = 0, byid = TRUE)),
            contour = FALSE, margin = FALSE, par.settings = rgrTheme, 
            maxpixels = pixel_num, colorkey = TRUE, main = "% Change in Relative Growth Rate") +
    levelplot(hill,  margin = FALSE, par.settings=hsTheme, maxpixels = pixel_num) +
    latticeExtra::layer(sp.polygons(cali_outline, lwd=2)) +
    latticeExtra::layer(sp.polygons(lobata_range,fill = "black", alpha = 0.15))
  
  print(p)
  
  # ### Plot of areas where SEs don't overlap
  # p = levelplot(rgr_rast_SEoverlap, 
  #               contour = FALSE, margin = FALSE,
  #               maxpixels = pixel_num, colorkey = TRUE, at = seq(0, 2, 1),
  #               main = "% Change in Relative Growth Rate") +
  #   levelplot(hill,  margin = FALSE, par.settings=hsTheme, maxpixels = pixel_num) +
  #   latticeExtra::layer(sp.polygons(cali_outline, lwd=2)) +
  #   latticeExtra::layer(sp.polygons(lobata_range,fill = "black", alpha = 0.15))
  # 
  # print(p)
  
  
  
  ## Histogram of changes in relative rgr
  
  ## Across all of california
  rgr_change_range = unlist(raster::extract(rgr_rast_change, lobata_range))
  rgr_change_cali = unlist(raster::extract(rgr_rast_change, 1:ncell(rgr_rast)))
  
  cat("Summary of rgr change within Valley oak range \n")
  print(summary(rgr_change_range))
  
  cat("Summary of rgr change across California \n")
  print(summary(rgr_change_cali))
  
  ## Growth within range
  p = ggplot(data.frame(rgr_change_range),
         aes(rgr_change_range)) + 
    geom_histogram(color = "black", fill = "steelblue2") +
    geom_vline(xintercept = 0, lty = 2) + 
    xlab("% change in relative growth rate") + 
    ggtitle("Within Valley oak range") +
    theme_bw() +
    theme(axis.text=element_text(size=12),
                axis.title=element_text(size=14),
                plot.margin = unit(c(1,1,1,1), "cm"))
  
  print(p)
  
  ## Growth within Cali
  p = ggplot(data.frame(rgr_change_cali),
         aes(rgr_change_cali)) + 
    geom_histogram(color = "black", fill = "steelblue2") + 
    geom_vline(xintercept = 0, lty = 2) +
    xlab("% change in relative growth rate") + 
    ggtitle("Across California") +
    theme_bw() +
    theme(axis.text=element_text(size=12),
                axis.title=element_text(size=14),
                plot.margin = unit(c(1,1,1,1), "cm"))
  
  print(p)
  
} # End function
  
  
  
  
  ## CMD
  plot_future_growth(fit_CMD$gam)
  
  ## Tmax
  plot_future_growth(fit_Tmax$gam)
  plot_future_growth(fit_Tmax_elev$gam)
  
  ## Tmin
  plot_future_growth(fit_Tmin$gam)
  
  ## DD5
  plot_future_growth(fit_DD5$gam)
  
  ## PCA model
  plot_future_growth(fit_clim_pca$gam, PCA = TRUE)
  
  

## HEIGHT PREDICTIONS - Use PC scores to make predictions about height




