#' Bisect Line Segments
#' 
#' Bisects each line segment in a polyline
#' @export
#' @param rpolyline depth values along a transect
#' @param projcrs
#' @param w width (m)
#' @param resolution resolution of bathymetry data (m)
#' @details Calculates the angle of the river centerline and draws a perpedicular bisector passing through it.
#' bisect_line_segments2 is preferred for the latest version of auto_transects
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


# ------------------------------------------------------------------------------------------------------------

#' Bisect line segments 2
#' 
#' Improved version of bisect_line_segments that uses the polyline_splitter from jmt2080ad's github page
#' @param rivsplit river centerline, split into segments
#' @param w width to extend the bisectors on either side of the river centerline
#' @param resolution resolution of bathymetry data (m)
#' @param mid TRUE to use the middle line segment, FALSE to assume the line segment is straight
#' @details When using mid = FALSE, this function implicitly assumes the line segments are not curved. 
#' The consequence of this is that the bisectors may not be somewhat off center, 
#' especially when the segments are relatively long. 
#' If line segments are relatively long, it is better to set mid = TRUE.
#' 
#' @example cross_section <- bisect_line_segments2(rivsplit, w = 1000, resolution = 5, mid = TRUE)
#' @export

bisect_line_segments2 <- function(rivsplit, w, resolution, mid = FALSE)
{
  
  # Initialize
  w <- w/resolution # width in pixels, resolution is 5 m, numerator is in meters 
  nseg <- length(rivsplit)
  lines_obj <- vector("list", length = nseg)
  coords.df <- array(dim = c(nseg, 4))
  coords.df <- as.data.frame(coords.df)
  colnames(coords.df) <- c("x1","y1","x2","y2")
  projcrs <- crs(rivsplit)
  
  for(seg in 1:nseg)
  {
    
    # Find the approximate midpoint
    # midpoints <- maptools::SpatialLinesMidPoints(rivsplit)
    
    nl <- length(rivsplit@lines[[seg]]@Lines)
    
    if (mid)
    {
      mid.ind <- round(nl/2)
      if (mid.ind == 0) # error catch
      {
        mid.ind <- 1
      }
      cm <- rivsplit@lines[[seg]]@Lines[[mid.ind]]@coords # coordinates of the middle line segment
      p1 <- data.frame(x = cm[1,1], y = cm[1,2])
      p2 <- data.frame(x = cm[2,1], y = cm[2,2])
    } else
    {
      c1 <- rivsplit@lines[[seg]]@Lines[[1]]@coords # coordinates of the first line segment
      c2 <- rivsplit@lines[[seg]]@Lines[[nl]]@coords # coordinates of the last line segment
      p1 <- data.frame(x = c1[1,1], y = c1[1,2]) # Get the first and last point
      p2 <- data.frame(x = c2[2,1], y = c2[2,2])
    }
  
    # Get the angle
    theta <- atan((p2$y - p1$y)/(p2$x - p1$x))
    if(is.na(theta)) # in case of vertical segment
    {
      theta <- pi/2
    } 
    theta_perp <- theta + pi/2
    
    # Draw the bisector
    xb1 <- (p1$x + p2$x)/2
    yb1 <- (p1$y + p2$y)/2
    
    # Extend it out across the river channel
    # (just needs to be at least half as wide as the widest part of the river)
    xb0 <- xb1 - w*resolution*cos(theta_perp)
    xb2 <- xb1 + w*resolution*cos(theta_perp)
    yb0 <- yb1 - w*resolution*sin(theta_perp)
    yb2 <- yb1 + w*resolution*sin(theta_perp)
    
    coords <- cbind(c(xb0, xb2), c(yb0, yb2)) # each line is defined by two points (OK)
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
