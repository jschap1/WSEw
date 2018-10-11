#' Estimate Widths
#' 
#' Estimates bankfull widths of transects based on the angle of the transect and the number of pixels in the transect.
#' @export
#' @param rpolyline a polyline crossing the river channel
#' @param tlength 
#' @param nseg 
#' @param resolution resolution of bathymetry data (m)
#' @details #' A simplification. Hopefully, it is accurate over large enough transects. 
#' Should probably revisit and make sure this code is reasonably accurate.
#' @return cross_section_width approximate bankfull width of the cross section
#' @examples 
#' wbf <- estimate_widths(rpolyline, resolution = 5, channel.pix, nseg)
#' @keywords bankfull width
#' @import sp

estimate_widths <- function(rpolyline, tlength, nseg, resolution = 5)
{
  cross_section_width <- vector(length=nseg)
  for (seg in 1:nseg)
  {
    cross_section_width[seg] <- tlength[[seg]]*resolution
  }
  return(cross_section_width)
}

# --------------------------------------------------------------------------------------------------------

# Get the first and last points on the line overlapping the main channel
# Make it a function, and do a quick analysis to see if it works.

#' Extract bankfull width
#' 
#' Uses first and last overlapped pixels to estimate bankfull width from 
#' perpendicular bisectors to the river centerline.
#' @param cross_section output from bisect_line_segments()
#' @param depth bathymetry raster
#' @examples wbf <- extract_wbf(cross_section[1:3], depth)
#' wbf_extract_wbf <- extract_wbf(cross_section, depth)
#' @export
#' @details It is slower than estimate_widths, but it is considerably more accurate.
#' If any of the transects are empty, the function will fail.

extract_wbf <- function(cross_section, depth)
{
  n.xs <- length(cross_section)
  wbf <- vector(length = n.xs)
  for (x in 1:n.xs)
  {
    transect <- extract(depth, cross_section[x], along = TRUE, cellnumbers = TRUE)
    
    main_channel <- get_main_channel(transect[[1]][,2], return_index = TRUE)
  
    first_cell <- transect[[1]][main_channel$first_ind,1]
    last_cell <- transect[[1]][main_channel$last_ind,1]
    
    p1 <- xyFromCell(depth, first_cell)
    p2 <- xyFromCell(depth, last_cell)
    
    wbf[x] <- dist(rbind(p1,p2)) 
    # estimate of bankfull width, ought to be subject to at most two 
    # pixels worth of error (though it can actually have more error than that)
    
    print(paste("Estimated wbf for cross section", x, "of", n.xs))
  }
  return(wbf)
}

# ------------------------------------------------------------------------------------------------------------

#' Extract bankfull width along transect
#' 
#' Method for calculating bankfull width using the package inlmisc from the USGS
#' @param cross_section output from bisect_line_segments()
#' @param depth bathymetry raster
#' @examples wbf1 <- extract_wbf_along_transect(cross_section[1:3], depth)
#' wbf_along_transect <- extract_wbf_along_transect(cross_section[1:42], depth)
#' @details Also fairly slow because it runs ExtractAlongTransect over each cross section, 
#' but it is inaccurate compared to extract_wbf
#' @export
#' @importFrom inlmisc

extract_wbf_along_transect <- function(cross_section, depth)
{
  
  n.xs <- length(cross_section)
  wbf <- vector(length = n.xs)
  
  for (x in 1:n.xs)
  {
    transect <- ExtractAlongTransect(cross_section[x], r = depth)
    main_channel <- get_main_channel(t1[[1]]$p21_depth, return_index = TRUE)
    
    if(all(is.na(main_channel))) # error catch in case transect is empty
    {
      wbf[x] <- NA
      next
    }
    
    wbf[x] <- transect[[1]]$dist[main_channel$last_ind] - transect[[1]]$dist[main_channel$first_ind]
    print(paste("Estimated wbf for cross section", x, "of", n.xs))
  }
  
  return(wbf)
}

# ------------------------------------------------------------------------------------------
# Compare the three methods

# measured_wbf <- vector(length = length(wbf_extract_wbf))
# measured_wbf[1:6] <- c(725,601,546,683,702,632)
# measured_wbf[37] <- 920
# measured_wbf[40:42] <- c(1098,865,946)
# 
# plot(wbf_estimate_widths, type = "l", main = "Comparing bankfull width estimation methods",
#      xlab = "cross section", ylab = "wbf (m)")
# points(wbf_estimate_widths, pch = 19, cex = 0.5)
# lines(wbf_extract_wbf, col = "red")
# points(wbf_extract_wbf, col = 'red', pch = 19, cex = 0.5)
# lines(wbf_along_transect, col = 'blue')
# points(wbf_along_transect, col = 'blue', pch = 19, cex = 0.5)
# points(measured_wbf, col = 'darkgreen', pch = 19, cex = 2)
# legend("topleft", 
#        legend = c("estimate_widths", "extract_wbf", "wbf_along_transect", "measured_wbf"), 
#        col = c("black","red","blue", "darkgreen"), lwd = c(1,1,1))
# 
# # save perpendicular bisectors as a shapefile
# # must convert to SpatialPolygonsDataFrame
# cross_section_sldf <- SpatialLinesDataFrame(cross_section, data = data.frame(a = 1:46))
# writeOGR(cross_section_sldf, 
#          dsn = "./Outputs/Cross_Sections/p21_transects", 
#          layer = "p21_transects", 
#          driver = "ESRI Shapefile")

# wbf_extract is (by far) the best method for estimating widths, based on comparison with manual measurements
# there is still a problem, though. The bisectors are not necessarily truly perpendicular to the river flowline
# also, they are not really evenly spaced, though this is less of a problem

# ------------------------------------------------------------------------------------------
# Scrap

# cross_section[1]@lines
# 
# depth[first_cell]

# Get xy location of the first and last cells
# Calculate linear distance between them
# Compare to width calculated with estimate_width and with USGS::ExtractAlongTransect
# This method is relatively slow.

# get index of t1 corresponding to beginning and end of the cross section
# write modified version of get_main_channel

# estimate_widths <- function(rpolyline, tlength, nseg, resolution = 5)
# {
#   cross_section_width <- vector(length=nseg)
#   for (seg in 1:nseg)
#   {
#     theta <- atan((rpolyline$y[seg+1] - rpolyline$y[seg])/(rpolyline$x[seg+1] - rpolyline$x[seg]))
#     if(is.na(theta)) {theta <- pi/2} # in case of vertical segment
#     theta_perp <- theta + pi/2
#     cross_section_width[seg] <- tlength[[seg]]*abs(cos(theta_perp))*resolution
#   }
#   return(cross_section_width)
# }

# if(all(is.na(main_channel$main_channel))) # error catch in case transect is empty
# {
#   wbf[x] <- NA
#   next
# }

# r <- 10
# transects1 <- extract(depth, cross_section[r], progress = "text", along = TRUE)
# transects2 <- extract(depth, cross_section[r], progress = "text")
# 
# par(mfrow =c(1,2))
# 
# plot(transects1[[1]], main = "along = true", type="l")
# plot(transects2[[1]], main = "along = false", type="l")
