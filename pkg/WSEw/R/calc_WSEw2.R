#' Calculate WSE-w
#' 
#' Calculates WSE-w pairs for each cross section.

#' @export
#' @param cross_sections list of x coordinates, bed elevations, and depths for each transect
#' @param interval number of points to sample in each cross-section. Unlike calc_WSE, interval does not determine the number of data points, which are instead constrained by the SWOT resolution (50m) and the river width
#' @param dx used for computing width with trapezoidal rule
#' @details 
#' Starts at w = 100m, the minimum value, and increases until bankfull width in 50 m increments
#' Not used. Use calc_WSEw, instead.
#' @return WSEw data frame containing WSE and flow width values
#' @examples 
#' xWSEw <- calc_WSEw(cross_sections, interval = 0.05, dx = 1)

# OUTPUTS
# WSE-w pairs for each cross section: xWSEw, cross section geometry

calc_WSEw2 <- function(cross_sections, interval = 0.05, dx = 1)
{
  
  x <- cross_sections$x
  b <- cross_sections$b
  d <- cross_sections$d
  
  # Pre-processing
  nseg <- length(b)
  dbf <- unlist(lapply(d, max, na.rm=TRUE)) # bankfull depth (m)
  wbf <- unlist(lapply(x, max, na.rm=TRUE)) # bankfull width (m)
  b.min <- unlist(lapply(b, min)) # bottom depth
  
  xWSEw <- vector(length = nseg, "list")
  for (seg in 1:nseg)
  {
    
    # Fit a linear spline to the observed channel geometry
    f1 <- approxfun(x = x[[seg]], y = b[[seg]], method="linear",
                    yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)
    
    # Make a vector of WSE(t) values to enter in with
    if (dbf[seg] < 0){next} # error check
    WSE <- seq(b.min[seg], b.min[seg] + dbf[seg], by = interval*dbf[seg])
    
    # Estimate the flow width for each WSE value
    w <- get_width(WSE, f1, x.max = max(x[[seg]]), delx = dx)
    
    # Use the WSE-w pairs to get WSE-w for widths ranging from 100m to wbf in 50 m increments
    approx.vals <- approx(w, y = WSE, seq(100, wbf[seg], by = 50), method="linear")
    WSEw <- data.frame(WSE = approx.vals$y, w = approx.vals$x)

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
