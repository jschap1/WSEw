# Identifying the channel centerline from the narivs dataset
#
# June 28,2018

# ------------------------------------------------------------------------
# Load libraries and data

library(raster)
library(rgdal)
library(maptools)

# Load bathymetry data
bathy.dir <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC"
bathy.name <- "bath_la_grange/bath_1997_lag/w001001.adf"
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
# /bath_la_grange/bath_1997_lag

# Pre-processing (gdal)
# Ensure narivs is in the desired coordinates (epsg 4326)
# Crop narivs to the extent of the depth raster for a given pool
# See narivs_preprocessing.sh for commands.

# Load river centerlines
riv.dir <- "/Users/jschap/Documents/Data/HydroSHEDs/na_riv_15s/Pools/"
riv <- readOGR(file.path(riv.dir, "nariv_cropped_pLAG.shp"))
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
# 1, 7, 9, 13, 17, 19 # Pool 7
# 1, 4, 5, 7, 10, 12, 17, 19, 24, 26 # Pool 8
# 4, 15, 20, 21, 23, 26, 29, 32, 35, 37, 41  # Pool 9
# 1, 5, 7, 11, 13, 15, 19, 21, 23, 25, 28, 29, 31, 34, 36 # Pool 10
# 1, 3, 5, 9, 10, 13, 16, 17, 20, 28, 32, 39, 57, 61, 64, 69 # Pool 13
# 15, 17, 20, 26, 27, 29, 30, 32 # Pool 26
#  # La Grange
# 8, 26, 29, 30, 31, 50, 59, 66, 71, 80, 85, 88, 100, 111, 115, 135, 141, 149, 152, 158, 166, 169, 170, 176, 188, 205, 215, 226, 240 

# Plot the main channel, segment by segment
linenums <- c(8, 31, 30, 26, 29, 50, 59, 66, 71, 80, 85, 88, 100, 111, 115, 135, 141, 149, 152, 158, 166, 170, 169, 176, 188, 205, 215, 226, 240)

plot(depth)
for (l in linenums)
{
  lines(riv@lines[[l]]@Lines[[1]]@coords)
}

i <- 20
i <- i+1
lines(riv@lines[[linenums[i]]]@Lines[[1]]@coords)

# 166, 170, 169, 176
# 8, 31, 30, 26, 29, 50
# get rid of 
# 26, 29, 169, 

# Make a list of Line objects using a loop or lapply
riv.lines.list <- vector(length = length(linenums), "list")
ind <- 1
for (r in linenums)
{
  riv.lines.list[[ind]] <- Lines(riv@lines[[r]]@Lines[[1]], ID = r) # Line
  ind <- ind + 1
}

# Put them all together
# centerline_p8 <- SpatialLines(riv.lines.list, proj4string = projcrs)

# Save the centerline
# save(centerline_p7, file = "Centerlines/centerline7.rda")

# --------------------------------------------------------------------------------------------------------

# Need to consolidate to have one feature only
coords <- array(dim = c(0,2))
for (l in linenums)
{
  coords <- rbind(coords, riv@lines[[l]]@Lines[[1]]@coords)
}

plot(coords)
# coords <- coords[-(10:30),]

# For Pool 26, need to re-order the coordinates from upstream to downstream
coords[,1] # easting
sortind <- order(coords[,1])
coords <- coords[sortind,]

# Make Spatial Lines object
coords <- as.data.frame(coords)
names(coords) <- c("x","y")
coordinates(coords) <- c("x","y")

# plot(coords)
# i <- 50
# i <- i + 1
# points(coords[i,], col = "red")

# Make Spatial Lines object
centerline_p26 <- SpatialLines(list(Lines(list(Line(coords)), "id")), proj4string = projcrs)

plot(depth)
# lines(centerline_p26)
lines(centerline_p26)

# Save the centerline (better format)
save(centerline_p26, file = "/Users/jschap/Desktop/Cross_Sections/Data/Centerlines/centerlineLAG.rda")
saveRDS(centerline_p26, file = "/Users/jschap/Desktop/Cross_Sections/Data/Centerlines/centerlineLAG.rds")
