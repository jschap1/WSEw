# Script to pre-process all the UMESC bathymetry data
# 10/5/2018

library(tools)
library(raster)
library(WSEw)

bathy.dir <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC"
setwd(bathy.dir)

# Bathymetry data file names:
bname <- vector(length = 9)
bname[1] <- file.path(bathy.dir, "/bath_pool_4/bath_1992_p4/w001001.adf")
bname[2] <- file.path(bathy.dir, "/pool_5/bath_2001_p5/w001001.adf")
bname[3] <- file.path(bathy.dir, "/bath_pool_7/bath_1998_p7/w001001.adf")
bname[4] <- file.path(bathy.dir, "/bath_pool_8/bath_1998_p8/w001001.adf")
bname[5] <- file.path(bathy.dir, "/bath_pool_9/bath_1999_p9/w001001.adf")
bname[6] <- file.path(bathy.dir, "/pool_10/bath_2001_p10/w001001.adf")
bname[7] <- file.path(bathy.dir, "bath_pool_13/bath_1997_p13/w001001.adf")
bname[8] <- file.path(bathy.dir, "/bath_pool_21/bath_1999_p21/w001001.adf")
bname[9] <- file.path(bathy.dir, "//bath_pool_26/bath_1997_p26/w001001.adf")

key <- c(4,5,7,8,9,10,13,21,26) # corresponding pool numbers

# Get raster depth values
for (i in 2:length(bname))
{
  
  if (i==8) next # skip pool 21, because those calculations have already been done
  umesc <- raster(bname[i])
  depth_5 <- as_numeric_raster(umesc)
  depth_5_refined <- depth_5
  depth_5_refined[values(depth_5>150)] <- NA
  depth_5 <- depth_5_refined
  writeRaster(depth_5, file = paste0("p", key[i], "_depth.tif"), progress = "text")
  print(paste("Saved", paste0("p", key[i], "_depth.tif")))
  
}

# Get cross section geometry

riv.dir <- "/Users/jschap/Desktop/Cross_Sections/Data/Centerlines"

rivname <- vector(length = 9)
rivname[1] <- file.path(riv.dir, "centerline4.rds")
rivname[2] <- file.path(riv.dir, "centerline5.rds")
rivname[3] <- file.path(riv.dir, "centerline7.rds")
rivname[4] <- file.path(riv.dir, "centerline8.rds")
rivname[5] <- file.path(riv.dir, "centerline9.rds")
rivname[6] <- file.path(riv.dir, "centerline10.rds")
rivname[7] <- file.path(riv.dir, "centerline13.rds")
rivname[8] <- file.path(riv.dir, "centerline21.rds")
rivname[9] <- file.path(riv.dir, "centerline26.rds")

# Reformatting centerline files so they are easier to loop through.
# i <- 8
# load(rivname[i])
# centerline <- centerline_p21
# saveRDS(centerline, file = file.path(riv.dir, "centerline21.rds"))

refWSE <- c(667, 660, 638.5, 631, 620, 611, 583, 470, 418.5)
halfwidth <- rep(1000, length(bname))
halfwidth[1] <- 3000
halfwidth[6] <- 2000
halfwidth[7] <- 2000

for (i in 1:length(bname))
{
  
  if (i==8) {next} # skip pool 21
  
  depth_5 <- raster(bname[i])
  transects_name <- paste0("p", key[i], "_xs_geometry.rda")
  riv <- readRDS(rivname[i])
  
  # Compute cross section data from raw bathymetry (transects_name is defunct)
  
  print(paste("Processing Pool", bname[i]))
  cross_sections <- auto_transects(section_length = 10e3, depth = depth_5, refWSE = refWSE[i],
                                   savename = transects_name, makeplot = TRUE, riv = riv, 
                                   halfwidth = halfwidth[i])
  
#  plot(depth_5)
#  lines(riv)
  
}