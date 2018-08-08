# WSEw project - code testing
# 
# Created 8/7/2018 JRS

# --------------------------------------------------------------------------------------------------
# Set up environment and load data

rm(list=ls())
#library(raster)
library(WSEw)
setwd("/Users/jschap/Desktop/Cross_Sections")

bathy.dir <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC"
bathy.name <- "bath_pool_21/bath_1999_p21"
bathyfile <- file.path(bathy.dir, bathy.name)

umesc <- raster(bathyfile)
levels(umesc)
depth_5 <- raster("Data/p21_depth.tif")

refWSE <- 470 # refWSE for pool 21 (ft)
refWSE <- refWSE/(39.37/12) #  convert to m

riv.dir <- "/Users/jschap/Desktop/Cross_Sections/Data/Centerlines"
load(file.path(riv.dir, "centerline21.rda"))
riv <- centerline_p21

saveloc <- "/Users/jschap/Desktop/Cross_Sections/Data/Processed_Data"
load(file.path(saveloc, "processed_data_p21_sl_5m_hires.rda"))
save(cross.sections, file = file.path(saveloc, "cross_sections_format2.rda"))

# -------------------------------------------------------------------------------------------------
# calc_WSEw / get_width

#cross.section <- list(x = cross_sections$x[[1]], b = cross_sections$b[[1]], d = cross_sections$d[[1]]) 
# may want to update to a more useful format, instead of a list of length 3, better to be a list of length nseg,
# except this would interfere with other functions that are set up for cross_sections to be in a particular format

# Using the updated format (format2) for the cross sections:
#x <- cross.sections[[1]]$x
#b <- cross.sections[[1]]$b
#d <- cross.sections[[1]]$d

#WSEw <- calc_WSEw(cross.section, interval = 0.05, dx = 1)




