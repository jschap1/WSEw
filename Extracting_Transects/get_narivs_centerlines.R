# Identifying the channel centerline from the narivs dataset
#
# June 28,2018

# ------------------------------------------------------------------------
# Load libraries and data

rm(list=ls())
library(raster)
library(rgdal)
library(maptools)

setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections")

# Load bathymetry data
bathy.dir <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC"
bathy.name <- "/pool_5/bath_2001_p5/w001001.adf"
depth <- raster(file.path(bathy.dir, bathy.name))
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

# Pre-processing (gdal)
# Ensure narivs is in the desired coordinates (epsg 4326)
# Crop narivs to the extent of the depth raster for a given pool
# See .txt file

# Load river centerlines
riv.dir <- "/Users/jschap/Documents/Data/HydroSHEDs/na_riv_15s/Pools/"
riv <- readOGR(file.path(riv.dir, "nariv_cropped_p5.shp"))
riv <- spTransform(riv, projcrs)

plot(depth, xlab = "Easting", ylab = "Northing", legend = FALSE)
lines(riv)

# --------------------------------------------------------------------------------------------------------
# Generate main channel centerline SpatialLines object

##### SKIP #####
linenum <- 23

# Get a line
riv@lines[[linenum]]@Lines

# Get coordinates of a line
riv@lines[[linenum]]@Lines[[1]]@coords

# Find longest line
ll <- vector(length = 48)
for (linenum in 1:48)
{
  ll[linenum] <- LineLength(riv@lines[[linenum]]@Lines[[1]]@coords)
}
################

# Issue: in narivs, there is not one continuous line for the main channel; they are divided up
# Soln: Go through one at a time and make a set of lines making up the centerline of the main channel
plot(depth)
linenum <- 1

linenum <- linenum+1
lines(riv@lines[[linenum]]@Lines[[1]]@coords)

# put the linenums of the main channel lines here:
# 5, 15, 19, 22, 23, 28, 34, 40, 45, 49, 52, 54, 59, 63, 67, 71, 73, 91 # Pool 4
# 2, 3, 6, 11, 14, 17 # Pool 5

# Plot the main channel, segment by segment
linenums <- c(2, 3, 6, 11, 14, 17)
plot(depth)
for (l in linenums)
{
  lines(riv@lines[[l]]@Lines[[1]]@coords)
}

# Make a list of Line objects using a loop or lapply
riv.lines.list <- vector(length = length(linenums), "list")
ind <- 1
for (r in linenums)
{
  riv.lines.list[[ind]] <- Lines(riv@lines[[r]]@Lines[[1]], ID = r) # Line
  ind <- ind + 1
}

# Put them all together
centerline_p5 <- SpatialLines(riv.lines.list, proj4string = projcrs)

# Save the centerline
save(centerline_p5, file = "Centerlines/centerline5.rda")

# --------------------------------------------------------------------------------------------------------

# Need to consolidate to have one feature only
coords <- array(dim = c(0,2))
for (l in linenums)
{
  coords <- rbind(coords, riv@lines[[l]]@Lines[[1]]@coords)
}

# Make Spatial Lines object
coords <- as.data.frame(coords)
names(coords) <- c("x","y")
coordinates(coords) <- c("x","y")

# Make Spatial Lines object
centerline_p5 <- SpatialLines(list(Lines(list(Line(coords)), "id")), proj4string = projcrs)

# Save the centerline (better format)
save(centerline_p5, file = "Centerlines/centerline5.rda")

