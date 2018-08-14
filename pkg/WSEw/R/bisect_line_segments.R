#' Bisect Line Segments
#' 
#' Bisects each line segment in a polyline
#' @export
#' @param rpolyline depth values along a transect
#' @param projcrs
#' @param w width (m)
#' @param resolution resolution of bathymetry data (m)
#' @details Calculates the angle of the river centerline and draws a perpedicular bisector passing through it.
#' @return cross_section a SpatialLines object containing all the bisectors
#' @examples 
#' bisect_line_segments(rpolyline, projcrs, w, resolution)
#' @keywords transects, hydraulics, cross sections
#' @import sp

bisect_line_segments <- function(rpolyline, projcrs, w, resolution)
{
  
  w <- w/resolution # width in pixels, resolution is 5 m, numerator is in meters 
  
  nseg <- dim(rpolyline)[1]-1
  coords.df <- array(dim = c(nseg, 4))
  coords.df <- as.data.frame(coords.df)
  colnames(coords.df) <- c("x1","y1","x2","y2")
  
  lines_obj <- vector("list", length = nseg)
  
  for (seg in 1:nseg)
  {
    theta <- atan((rpolyline$y[seg+1] - rpolyline$y[seg])/(rpolyline$x[seg+1] - rpolyline$x[seg]))
    if(is.na(theta)) {theta <- pi/2} # in case of vertical segment
    theta_perp <- theta + pi/2
    
    # Draw the bisector
    xb1 <- (rpolyline$x[seg] + rpolyline$x[seg+1])/2
    yb1 <- (rpolyline$y[seg] + rpolyline$y[seg+1])/2
    
    # Extend it out across the river channel
    # (just needs to be at least half as wide as the widest part of the river)
    xb0 <- xb1 - w*resolution*cos(theta_perp)
    xb2 <- xb1 + w*resolution*cos(theta_perp)
    yb0 <- yb1 - w*resolution*sin(theta_perp)
    yb2 <- yb1 + w*resolution*sin(theta_perp)
    
    coords <- cbind(c(xb0, xb2), c(yb0, yb2)) # each line is defined by two points
    coords.df[seg,] <- c(coords[1,1], coords[1,2], coords[2,1], coords[2,2])
    # Outputs a data frame with coordinates for the endpoints of the bisectors of rpolyline
    
    # Make a SpatialLines object for each bisector
    simple_line <- Line(cbind(rbind(coords.df$x1[seg], coords.df$x2[seg]),
                              rbind(coords.df$y1[seg], coords.df$y2[seg])))
    lines_obj[[seg]] <- Lines(list(simple_line), ID = as.character(seg)) # can list these
    
  }
  cross_section <- SpatialLines(lines_obj,proj4string = projcrs)
  return(cross_section)
}