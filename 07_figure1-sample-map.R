# Source files ------------------------------------------------------------

source("./01_clean-process-explore-data.R")
source("./03_adding-climate-data.R")


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

lobata_range <- readShapePoly("./data/gis/valley_oak_range/querloba.shp",
                              proj4string = CRS("+proj=longlat +datum=WGS84"))

all_moms <- SpatialPoints(climate_garden_mom[, c("longitude", "latitude")], 
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))

gbs_moms_sub <- climate_gbs_mom %>%
              dplyr::filter(climate_gbs_mom$accession %in% dat_all_scaled$accession)
plot(all_moms)

gbs_moms <- SpatialPoints(gbs_moms_sub[, c("longitude", "latitude")],
                          proj4string = CRS("+proj=longlat +datum=WGS84"))
plot(gbs_moms)

garden_sites <- SpatialPoints(garden_climate[, c("longitude", "latitude")],
                              proj4string = CRS("+proj=longlat +datum=WGS84"))
plot(garden_sites)




# Make and save Plot ---------------------------------------------------------------

     pixel_num = 1e5
      hsTheme <- modifyList(GrTheme(), list(regions=list(alpha=.15)))
       demTheme <- modifyList(RdBuTheme(), list(regions=list(alpha=.5)))
       
    p = levelplot(dem, contour = FALSE, margin = FALSE, par.settings = demTheme, 
                  maxpixels = pixel_num, colorkey = TRUE) +
      levelplot(hill,  margin = FALSE, par.settings=hsTheme, maxpixels = pixel_num) +
      latticeExtra::layer(sp.polygons(cali_outline, lwd=2)) + 
      latticeExtra::layer(sp.polygons(lobata_range,fill = "black", alpha = 0.25)) +
      latticeExtra::layer(sp.points(all_moms, pch = 21, cex = 1.5,
                                    fill = "steelblue2", col = "grey10", alpha = .5)) +
      latticeExtra::layer(sp.points(gbs_moms, pch = 21, cex = 1.5,
                                    fill = "forestgreen", col = "grey10", alpha = .5))  +
      latticeExtra::layer(sp.points(garden_sites, pch = 22, cex = 2, 
                                    col = "grey10", fill = "orange"))
    
    print(p)
    
    pdf(paste0("./figs_tables/Figure1_sample_map_", Sys.Date(), ".pdf"))
    print(p)
      dev.off()

