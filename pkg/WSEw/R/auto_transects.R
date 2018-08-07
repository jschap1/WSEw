#' Auto Transects
#' 
#' Builds cross sections from bathymetry data
#' @param section_length spatial discretization for dividing the river into cross sections (m)
#' @param riv polyline file for the main river channel
#' @param depth raster of depth values (m)
#' @param refWSE Reference WSE for the measured bathymetry data (m)
#' @param savename where to save the outputs
#' @param halfwidth guess for how wide the channel is. Default is 1000 m.
#' @param k smoothing width
#' @param makeplot whether or not to make a plot. Default is FALSE.
#' @keywords transects, hydraulics, cross sections
#' @import sp
#' @importFrom raster extract
#' @export
#' @examples 
#' auto_transects(section_length, savename = "p21.rda")

auto_transects <- function(section_length, riv, depth, refWSE, 
                           savename, halfwidth = 1000, k = 5, 
                           makeplot = FALSE)
{
  
  # Draw transects and extract depths
  # depth <- raster(bathy)
  # depth[depth==9999] <- NA # 9999 is a code that means something, but remove it for our purposes
  projcrs <- crs(depth)
  
  # Creates a single "polyline" file for input to resample_polyline
  x <- coordinates(riv@lines[[1]])[[1]][,1]
  y <- coordinates(riv@lines[[1]])[[1]][,2]
  polyline = data.frame(x = x, y = y)
  
  # Divides the polyline into equal-length segments
  rpolyline <- resample_polyline(polyline, interval_length = section_length)
  
  # Bisect line segments
  cross_section <- bisect_line_segments(rpolyline, projcrs, halfwidth, resolution = res(depth)[1])
  
  # Plot to check
  if (makeplot)
  {
    plot(depth, main = "UMRB Bathymetry, Pool 4", xlab = "Easting", ylab = "Northing", legend = FALSE)
    lines(riv)
    lines(cross_section, col = "Red") # the segments are numbered from north to south
  }

  print("Extracting values along transects, this can take a long time (10-60 minutes)")
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
  
  # d <- lapply(main_channel, get_depth_from_lutable, depth)
  d <- main_channel
  
  # Assume the first and last point are zero depth (banks)
  for (seg in 1:nseg)
  {
    d[[seg]][1] <- 0
    d[[seg]][channel.pix[seg]] <- 0
  }
  
  # cross section bankfull width (not entirely correct)
  wbf <- estimate_widths(rpolyline, resolution = res(depth), channel.pix, nseg) 
  
  x <- vector("list", length = nseg) # x coordinate, using river banks as beginning and end
  for (seg in 1:nseg)
  {
    x[[seg]] <- seq(0, wbf[seg], length.out = length(d[[seg]]))
  }
  
  b <- lapply(d, function(x, refWSE) {refWSE-x}, refWSE)
  
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
    d.smooth[[seg]] <- filter(d[[seg]], sides = 2, filter = rep(1/k,k))
  }
  
  # Remove empty channels
  if (any(channel.pix<=10))
  {
    rm.ind <- which(channel.pix<=10)
    d <- d[-rm.ind]
    b <- b[-rm.ind]
    x <- x[-rm.ind]
  }
  
  # Save extracted transect information
  save(x, b, b.smooth, d, d.smooth, channel.pix, file = savename)
  print(paste("Saved transect data to", savename))
  
  cross.sections <- list(x = x, b = b, d = d)
  return(cross.sections)
  
}
