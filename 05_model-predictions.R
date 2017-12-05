


##########################
### RASTER PROCESSING ####
##########################

## Read in elevation dem and create a hillshade
dem <- raster("../data/dem/ClimateNA_DEM_cropped.tif")

slope = terrain(dem, opt='slope')
aspect = terrain(dem, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)


## Get general file names without extension
general_files <- list.files(path="../data/climate data/NA_NORM_8110_Bioclim_netCDF/",
                            pattern = "_cropped.tif$", full.names = FALSE)
general_files <- gsub(pattern = ".tif", replacement = "", general_files)

var_names <- paste(gsub(pattern = "_cropped", "", general_files), "_dif", sep = "")

general_files <- paste("/", general_files, sep = "") ## Need this for grep to be specific

## Read full files names of current climates
current_files <- list.files(path="../data/climate data/NA_NORM_8110_Bioclim_netCDF/",
                            pattern = "_cropped.tif$", full.names = TRUE)

## Read full file names of future climate
future_files <- list.files(path="../data/climate data/NA_ENSEMBLE_rcp85_2080s_Bioclim_netCDF/",
                           pattern = "_cropped.tif$", full.names = TRUE)


## Loop through climate variables and calculate difference between future and current climate

## Initialize stack of rasters  
dif_rast <- stack()
cur_rast <- stack()
x = 1 # For looping through var names in plot title
for(var in general_files){
  
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
  
  ## Plot output
  # plot(difference_rast, main = paste("Difference in:", var_names[x]))
  
  ## Add to stack of rasters
  dif_rast <- raster::stack(dif_rast, difference_rast)
  cur_rast <- raster::stack(cur_rast, current_rast)
  x = x+1
}  

names(dif_rast) <- var_names
names(cur_rast) <- gsub("_dif", "", var_names)

dif_rast_raw  <- dif_rast




#dif_rast_raw <- stack(dif_rast, cur_rast)





## Cropping based on lobata range
dif_rast = dif_rast_raw
# dif_rast <- raster::mask(dif_rast_raw, rgeos::gBuffer(lobata_range,                                                      width = 0, byid = TRUE))


## Extract values from raster into dataframe
extracted_vals = extract(dif_rast, 1:ncell(dif_rast))
extracted_vals <- as.data.frame(extracted_vals)


### Scale extracted values to match with PCA

extracted_vals_scaled <- extracted_vals

for(var in climate_vars_dif){
  extracted_vals_scaled[, var] <- (extracted_vals[, var] - climate_var_means[var]) / 
    climate_var_sds[var]
}

summary(extracted_vals_scaled)

## Make histogram of predicted climate differences vs. those from the gardens
# 
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
# ## MCMT
# levelplot(dif_rast[["MCMT_dif"]], contour = FALSE, margin = FALSE, main = "MCMT_dif",
#           par.settings = viridisTheme()) +  layer(sp.polygons(cali_outline, lwd=2))  +
#   layer(sp.polygons(lobata_range, fill = "black", alpha = 0.35))
# 
# ## MWMT
# levelplot(dif_rast[["MWMT_dif"]], contour = FALSE, margin = FALSE, main = "MWMT_dif",
#           par.settings = viridisTheme()) +  layer(sp.polygons(cali_outline, lwd=2))  +
#   layer(sp.polygons(lobata_range, fill = "black", alpha = 0.35))
# 
# ## DD5
# levelplot(dif_rast[["DD5_dif"]], contour = FALSE, margin = FALSE, main = "DD5_dif",
#           par.settings = viridisTheme()) +  layer(sp.polygons(cali_outline, lwd=2))  +
#   layer(sp.polygons(lobata_range, fill = "black", alpha = 0.35))
# 


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
# 


#### Climate vars model - making predictions

predictions <- extracted_vals_scaled[, climate_vars_dif]

predictions$cluster_assigned = 1
predictions$height_2016_scaled = 0


## Extract and scale elevation data

predictions$elev <- extract(dem,  1:ncell(dem))
predictions$elev <- (predictions$elev - elev_mean) / elev_sd


## HEIGHT PREDICTIONS - Use PC scores to make predictions about height

height_predictions =  predict(fit_gen_climPC$gam, predictions, re.form = NA, type = "response")

summary(height_predictions)


height_rast <- dif_rast[[1]] ## Initialize raster

values(height_rast) <- c(height_predictions)


## Plot of predicted relative height

hsTheme <- modifyList(GrTheme(), list(regions=list(alpha=.15)))
heightTheme <- modifyList(RdBuTheme(), list(regions=list(alpha=1)))

pixel_num = 1e5 ## Can make resolution better by making this 1e6 

levelplot(height_rast, 
          #levelplot(raster::mask(height_rast, rgeos::gBuffer(lobata_range, width = 0, byid = TRUE)),
          contour = FALSE, margin = FALSE, par.settings = heightTheme, 
          maxpixels = pixel_num, colorkey = TRUE, main = "Relative height") +
  levelplot(hill,  margin = FALSE, par.settings=hsTheme, maxpixels = pixel_num) +
  layer(sp.polygons(cali_outline, lwd=2)) + layer(sp.polygons(lobata_range,
                                                              fill = "black", alpha = 0.15))


## Histogram of changes in relative height

## Across all of california

height_relative_range = unlist(raster::extract(height_rast, lobata_range))
height_relative_cali = unlist(raster::extract(height_rast, 1:ncell(height_rast)))

summary(height_relative_range)
summary(height_relative_cali)

## Growth within range
ggplot(data.frame(height_relative_range),
       aes(height_relative_range)) + 
  geom_histogram(color = "black", fill = "steelblue2") + theme_bw() +
  geom_vline(xintercept = 1, lty = 2) + ggtitle("Growth within range") + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

sum(height_relative_range > 1, na.rm = TRUE) / length(na.omit(height_relative_range))

## Growth within Cali
ggplot(data.frame(height_relative_cali),
       aes(height_relative_cali)) + 
  geom_histogram(color = "black", fill = "steelblue2") + theme_bw() +
  geom_vline(xintercept = 1, lty = 2) + ggtitle("Growth across California") + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

sum(height_relative_cali > 1, na.rm = TRUE) / length(na.omit(height_relative_cali))





### Calculating future temps within q lobata range

temp = raster("../data/climate data/NA_ENSEMBLE_rcp85_2080s_Bioclim_netCDF/Tave_sm_cropped.tif")

temp = raster("../data/climate data/NA_NORM_8110_Bioclim_netCDF/Tave_sm_cropped.tif")

temp = raster::mask(temp, rgeos::gBuffer(lobata_range,                                                      width = 0, byid = TRUE))

plot(temp)

summary(temp)


