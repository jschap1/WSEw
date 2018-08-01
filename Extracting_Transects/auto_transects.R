# UMRB transects
#
# Updated to use centerlines derived from the gridded bathymetry dataset 6/21/2018
# Cleaned up 6/22/2018
# Updated 6/26/2018

# ------------------------------------------------------------------------
# Load libraries, data, and helper functions

rm(list=ls())
library(raster)
library(rgdal)
library(maptools)

setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections")
source("resample_polyline.R")
source("auto_transects_helper_functions.R")

# Load bathymetry data
bathy.dir <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC"
bathy.name <- "bath_pool_21/bath_1999_p21"
depth <- raster(file.path(bathy.dir, bathy.name))
depth[depth==9999] <- NA # 9999 is a code that means something, but remove it for our purposes
projcrs <- crs(depth)

# Bathymetry data locations index:
# /bath_pool_4/bath_1992_p4
# /pool_5/bath_2001_p5
# /bath_pool_7/bath_1998_p7
# /bath_pool_8/bath_1998_p8
# /bath_pool_9/bath_1999_p9
# /pool_10/bath_2001_p10
# /bath_pool_13/bath_1997_p13
# /bath_pool_21/bath_1999_p21
# /bath_pool_26/bath_1997_p26

# Load river centerlines
riv.dir <- "Centerlines"
load(file.path(riv.dir, "centerline21.rda"))
riv <- centerline_p21

plot(depth, xlab = "Easting", ylab = "Northing", legend = FALSE)
lines(riv) # if I could extract the main channel from narivs, then I could just use that

# ------------------------------------------------------------------------
# Draw transects and extract depths

section_length <- 500 # m
n <- 1 # index for pool; this code is for only one pool at a time 

# Reference WSE
refWSE <- 470 # ft
refWSE <- refWSE*12/39.37

# Creates a single "polyline" file for input to resample_polyline
x <- coordinates(riv@lines[[n]])[[1]][,1]
y <- coordinates(riv@lines[[n]])[[1]][,2]
polyline = data.frame(x = x, y = y)

# Divides the polyline into equal-length segments
rpolyline <- resample_polyline(polyline, interval_length = section_length)

# Bisect line segments
w <- 1000 # meters, expected maximum half-width of channel, tune it
cross_section <- bisect_line_segments(rpolyline, projcrs, depth, w)

# Plot to check
plot(depth, main = "UMRB Bathymetry, Pool 4", xlab = "Easting", ylab = "Northing", legend = FALSE)
lines(riv)
lines(cross_section, col = "Red") # the segments are numbered from north to south

# after finding nice cross sections, can highlight them
lines(cross_section[nice.ind], col = "Blue") 
lines(cross_section[-nice.ind], col = "Red")
legend("topright", legend = c("Fitted","Not fitted"), col = c("Blue", "Red"), lty = c(1,1))

# would also like to label them, but unsure how
#spplot(cross_section)
#lineLabel(cross_section@lines, label=seq(1,103,by=1))

# Extract values along transects, this can take a long time (10-60 minutes)
transects <- extract(depth, cross_section, progress = "text")

# Remove null (empty) transects
null.ind <- unlist(lapply(transects, is.null))
na.ind <- lapply(transects, is.na)
na.ind <- unlist(lapply(na.ind, all))
if (sum(null.ind)>0)
{
  transects <- transects[-which(null.ind)]
}
if (sum(na.ind)>0)
{
  transects <- transects[-which(na.ind)]
}
nseg <- length(transects)

# ------------------------------------------------------------------------------
# Extract x-y information for plotting transects

# Get distance of each transect

main_channel <- lapply(transects, get_main_channel)

# Zero length channels cause problems
channel.pix <- unlist(lapply(main_channel, length)) # number of values in each transects list value
# not the same as channel width because of the angle

transects.depth <- lapply(main_channel, get_depth)

# Assume the first and last point are zero depth (banks)
for (seg in 1:nseg)
{
  transects.depth[[seg]][1] <- 0
  transects.depth[[seg]][channel.pix[seg]] <- 0
}

# tlength <- SpatialLinesLengths(main_channel, longlat=FALSE) 
# not valid because main_channel is shorter than transects
# Instead, make use of the fact that the data are 5 m resolution - need angle, too.
# Channel width is an important output, or is it? Using normalized widths.

# cross section bankfull widths (Fix this up)
widths <- estimate_widths(rpolyline, resolution = 5, channel.pix, nseg) 

tdist <- vector("list", length = nseg) # x coordinate, using river banks as beginning and end
for (seg in 1:nseg)
{
  tdist[[seg]] <- seq(0, widths[seg], length.out = length(transects.depth[[seg]]))
}

b <- lapply(transects.depth, function(x, refWSE) {refWSE-x}, refWSE)

# Smooth using a k-point moving average
k <- 5
b.smooth <- vector("list", length = nseg)
d.smooth <- vector("list", length = nseg)
for (seg in 1:nseg)
{
  if (channel.pix[seg] <= k)
  { # error handling for zero-width channels
    next
  }
  b.smooth[[seg]] <- filter(b[[seg]], sides = 2, filter = rep(1/k,k))
  d.smooth[[seg]] <- filter(transects.depth[[seg]], sides = 2, filter = rep(1/k,k))
}

# ------------------------------------------------------------------------------
# Save extracted transect information
save(tdist, b, b.smooth, transects.depth, d.smooth, widths, channel.pix, file ="Transects/pool_21_500m.rda")

# Also may want to save:
# northing, easting of each coordinate on the transect
# pool number, reference WSE
# chainage at each transect
# distance between transects

# Save bathymetry raster in a convenient format for sharing
depth.vals <- get_depth(values(depth)) # takes an hour or so
depth.copy <- depth
values(depth) <- depth.vals
writeRaster(depth, filename = "p21_depth_utm.tif")
# Convert to wgs84 in gdal.
