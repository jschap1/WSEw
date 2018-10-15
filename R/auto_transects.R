#' Auto Transects
#' 
#' Builds cross sections from bathymetry data
#' @export
#' @param section_length spatial discretization for dividing the river into cross sections (m)
#' @param riv polyline file for the main river channel
#' @param depth raster of depth values (m)
#' @param refWSE Reference WSE for the measured bathymetry data (m)
#' @param savename where to save the outputs
#' @param halfwidth guess for how wide the channel is. Default is 1000 m.
#' @param k smoothing width
#' @param makeplot whether or not to make a plot. Default is FALSE.
#' @details It is a pretty long function. Here is a synopsis:
#' 1. Breaks up a polyline of the river centerline into roughly evenly-spaced segments of a specified length
#' 2. Draws a perpendicular bisector to each segment. Choose a good value of halfwidth to ensure the bisector completely crosses the river.
#' 3. Optionally plots the river centerline, gridded bathymetry data, and the bisectors (henceforth called "cross-sections")
#' 4. Extracts depth values along the cross-sections using extract(). This is the most time-consuming step.
#' 5. Removes empty cross-sections, those for which there are no available bathymetry data
#' 6. Gets the main river channel in the case when there are multiple channels by selecting the widest channel
#' 7. Estimates bankfull widths at each cross section. The estimates may be inaccurate for narrow river segments.
#' 8. Uses reference WSE to calculate bed elevation along each cross section.
#' 9. Smooths bed elevations and depths using a k-point moving average, with k=5.
#' 10. Saves cross-section information as savename.
#' 11. Puts x, b, d (distance from horizontal datum, bed elevation, depth) in a list called cross_sections
#' @return cross_sections geometry data for cross sections along the river
#' @examples 
#' cross_sections <- auto_transects(section_length = 5, depth = depth, refWSE = refWSE, savename = transects_name, makeplot = FALSE, riv = riv)
#' @keywords transects, hydraulics, cross sections
#' @import sp
#' @importFrom raster extract

auto_transects <- function(section_length, riv, depth, refWSE, 
                           savename, halfwidth = 1000, k = 5, 
                           makeplot = FALSE)
{
  
  rivsplit <- splitLines(riv, dist = section_length)
  projcrs <- crs(depth)
  crs(rivsplit) <- projcrs
  
  # Write out the shapefile
  # riv.df <- SpatialLinesDataFrame(riv, data = data.frame(ID = 1), match.ID = FALSE)
  # writeOGR(riv.df, dsn = "riv", layer = "riv", driver = "ESRI Shapefile")
  
  # Bisect line segments
  cross_section <- bisect_line_segments2(rivsplit, w = halfwidth, resolution = res(depth)[1], mid = TRUE)
  
  # Plot to check
  if (makeplot)
  {
    plot(depth, main = "Bathymetry", xlab = "Easting", ylab = "Northing", legend = FALSE)
    lines(riv.smooth.ksmooth)
    # lines(riv)
    lines(cross_section[40:50], col = "Red")
    lines(cross_section, col = "Red") # the segments are numbered from north to south
  }

  print("Extracting values along transects.")
  print("This can take a LONG time.")
  print("It takes about 12 seconds per cross section per processor")
  res <- extract_xs_wbf(cross_section[40:50], depth, hpc = FALSE)
  wbf <- res$wbf
  main_channel <- res$main_channel
  
  # Remove null (empty) transects
  null.ind <- unlist(lapply(main_channel, is.null))
  na.ind <- lapply(main_channel, is.na)
  na.ind <- unlist(lapply(na.ind, all))
  if (sum(null.ind)>0)
  {
    main_channel <- main_channel[-which(null.ind)]
  }
  if (sum(na.ind)>0)
  {
    main_channel <- main_channel[-which(na.ind)]
  }
  nseg <- length(main_channel)
  
  if(any(is.na(wbf)))
  {
    wbf.na.ind <- which(is.na(wbf))
    wbf <- wbf[-wbf.na.ind]
  }

  # ------------------------------------------------------------------------------
  # Extract x-y information for plotting transects
  
  # Get distance of each transect
  # main_channel <- lapply(transects, get_main_channel)
  
  # Zero length channels cause problems
  channel.pix <- unlist(lapply(main_channel, length)) # number of values in each transects list value
  # not the same as channel width because of the angle
  
  # d <- lapply(main_channel, get_depth_from_lutable, depth)
  d <- main_channel # depths
  
  # Assume the first and last point are zero depth (banks)
  for (seg in 1:nseg)
  {
    d[[seg]][1] <- 0
    d[[seg]][channel.pix[seg]] <- 0
  }
  
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
  
  # Remove (nearly) empty channels
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
  
  cross.sections <- list(x = x, b = b, d = d) # old format
  
  # this may be a more useful format
  #cross.sections <- vector(length = nseg, "list")
  #for (seg in 1:nseg)
  #{
  #  cross.sections[[seg]] <- data.frame(x = x[[seg]], b = b[[seg]], d = d[[seg]])
  #}
  
  return(cross.sections)
  
}
