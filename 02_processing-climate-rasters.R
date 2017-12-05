

## This script loads a series of raster files in a directory, crops them to california outline (more or less), projects them into lat / long, and saves them as a geotiff format

## Load required packages
library(raster)
library(sp)
library(maptools)


## Load in california outline
  cali_outline <- readShapePoly("./data/gis/california_outline/california_outline.shp",
                                      proj4string = CRS("+proj=longlat +datum=WGS84"))

## project to lamber conformal conic to use for cropping below - transforming this and then cropping speeds up operations a lot rather than transforming raster and cropping that

## CRS comes from read in .nc file with raster() and copying the auto-read CRS, adjusting for a typo
  cali_outline_lcc <- spTransform(cali_outline, CRS("+proj=lcc +lat_1=49.0 +lat_2=77.0 +lat_0=0.0 +lon_0=-95.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))


## Read in climate variables - https://adaptwest.databasin.org/pages/adaptwest-climatena

  dir_name <- ".data/gis/dem/"
  
  raster_files <- list.files(dir_name, full.names = TRUE)
  raster_files <- raster_files[grep("*.nc", raster_files)] # Only those with nc extension


## Loop over all rasters in the directory
  for(file in raster_files){
    
    cat("Working on file:", file, "... \n")
    
    ## Load in raster file
    raster_temp <- raster(file)
    
    ## Assign correct CRS since there is a typo in the NetCDF file
    crs(raster_temp) <- "+proj=lcc +lat_1=49.0 +lat_2=77.0 +lat_0=0.0 +lon_0=-95.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    
    ## Crop
    raster_temp_cropped <- raster::crop(raster_temp, cali_outline_lcc)
    
    ## Masking to get exact shape
    raster_temp_cropped <- mask(raster_temp_cropped, cali_outline_lcc)
    
    ## Project to latlong
    raster_cropped_latlong <- projectRaster(raster_temp_cropped, crs="+proj=longlat +datum=WGS84")
    
    ## Write to file - geotiff format set by setting file extension
    
    file_new = gsub(".nc", "_cropped.tif", file)
    
    writeRaster(raster_cropped_latlong, file_new, overwrite = TRUE)
    
  }

