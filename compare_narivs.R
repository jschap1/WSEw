#' Compare Andreadis et al. (2013) hydraulic geometry-based bathymetry dataset to UMESC data
#' 
#' 

library(rgdal)

nariv.name <- "/Users/jschap/Documents/Data/HydroSHEDs/na_riv_15s/Pools/na_riv_15s_wgs_cropped_p21.shp"

nariv <- readOGR(nariv.name)
