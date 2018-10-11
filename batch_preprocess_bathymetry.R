# Script to pre-process all the UMESC bathymetry data
# 10/5/2018

library(tools)
library(raster)
library(WSEw)

setwd("/Users/jschap/Documents/Research/SWOTBATH")
bathy.dir <- "./Data/UMESC_From_Box/"

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
bname[9] <- file.path(bathy.dir, "/bath_pool_26/bath_1997_p26/w001001.adf")

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

riv.dir <- "./Data/HydroSHEDs/Centerlines"

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

refWSE <- c(667, 660, 638.5, 631, 620, 611, 583, 470, 418.5)
halfwidth <- rep(1000, length(bname))
halfwidth[1] <- 3000
halfwidth[6] <- 2000
halfwidth[7] <- 2000

dname <- vector(length = 9) # names of depth rasters (actually is is bed elevation, but whatever)
for (i in 1:9)
{
  dname[i] <- file.path(bathy.dir, paste0("p", key[i], "_depth.tif"))
}

for (i in 1:length(bname))
{
  
  if (i==8) {next} # skip pool 21
  
  depth_5 <- raster(dname[i])
  transects_name <- paste0("./Outputs/Cross_Sections/p", key[i], "_xs_geometry.rda")
  transects_name <- "./Outputs/Cross_Sections/p21_test_smooth.rda"
  riv <- readRDS(rivname[i])
  
  # Smooth the river centerline
  # riv.smooth.chaikin <- smooth(riv, method = "chaikin")
  riv.smooth.ksmooth <- smooth(riv, method = "ksmooth", smoothness = 1)
  # riv.smooth.ksmooth <- smooth_ksmooth(riv, smoothness = 5)
  
  plot(depth_5)
  lines(riv)
  lines(riv.smooth.ksmooth, col = "blue")
  # plot(riv.smooth.ksmooth, col = "blue")
  
  # lines(riv.smooth.ksmooth, col = "red")
  
  # riv.smooth.spline <- smooth(riv, method = "spline")
  # plot(riv)
  # lines(riv.smooth.chaikin, col = "red")
  # lines(riv.smooth.ksmooth, col = "blue")
  # lines(riv.smooth.spline, col = "green")
  
  # Compute cross section data from raw bathymetry (transects_name is defunct)
  
  print(paste("Processing Pool", bname[i]))
  cross_sections <- auto_transects(section_length = 5e3, depth = depth_5, refWSE = refWSE[i],
                                   savename = transects_name, makeplot = TRUE, riv = riv.smooth.ksmooth, 
                                   halfwidth = halfwidth[i])
  
  xsname <- paste0("cross_sections_p", key[i], ".rds")
  saveRDS(cross_sections, file = file.path("./Outputs/Cross_Sections", xsname))
  
  # plot(depth_5)
  # lines(riv)
  r <- 2
  plot(cross_sections$x[[r]], cross_sections$b[[r]], type = "l", xlab="x", ylab = "b")
  
}