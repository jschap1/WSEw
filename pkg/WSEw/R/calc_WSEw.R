#' Calculate WSE-w
#' 
#' Calculates WSE-w pairs for each cross section
#' @param cross_sections list of x coordinates, bed elevations, and depths for each transect
#' @param interval number of points to sample in each cross-section
#' @param dx used for computing width with trapezoidal rule
#' @export
#' @examples calc_WSEw(x, b, interval = 0.05, dx = 1)

# OUTPUTS
# WSE-w pairs for each cross section: xWSEw, cross section geometry

calc_WSEw <- function(cross_sections, interval = 0.05, dx = 1)
{

  x <- cross_sections$x
  b <- cross_sections$b
  d <- cross_sections$d
  
  # Pre-processing
  nseg <- length(b)
  dbf <- unlist(lapply(d, max, na.rm=TRUE)) # bankfull depth (m)
  b.min <- unlist(lapply(b, min)) # bottom depth
  
  xWSEw <- vector(length = nseg, "list")
  for (seg in 1:nseg)
  {

    # Make a vector of WSE(t) values to enter in with
    if (dbf[seg] < 0){next} # error check
    WSE <- seq(b.min[seg], b.min[seg] + dbf[seg], by = interval*dbf[seg])
    
    # Estimate the flow width for each WSE value
    w <- get_width(WSE, x[[seg]], b[[seg]], delx = dx)
    
    # Put WSE-w info into a data frame
    WSEw <- data.frame(WSE = WSE, w = w)

    xWSEw[[seg]] <- WSEw
  }
  
  # Remove null cross sections, those with no measurements (optional)
  rm.ind <- which(unlist(lapply(xWSEw, is.null)))
  if (any(rm.ind))
  {
    xWSEw <- xWSEw[-rm.ind] # this doesn't always seem to work...
  }
  
  return(xWSEw)
}
