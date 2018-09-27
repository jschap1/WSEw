#' Generate h-w
#' 
#' Generates height and width data for each cross section, as it might be sampled by SWOT
#' @export
#' @param cross_sections list of x coordinates, bed elevations, and depths for each transect
#' @param dist type of distribution for SWOT-sampled stage values.
#' @param pars vector of distribution parameters
#' @param n.obs number of SWOT observations to make. The UMRB is observed 1-2 times per 21-day cycle
#' @param dx used for computing width with trapezoidal rule
#' @details Makes a vector of WSE values from the minimum bed elevation to the bankfull WSE
#' Computes flow width for each WSE-w value.
#' Puts the estimated flow widths into the WSEw data frame
#' @return WSEw data frame containing WSE and flow width values
#' @examples 
#' xWSEw <- calc_WSEw()

# OUTPUTS
# WSE-w pairs for each cross section: xWSEw, cross section geometry

generate_hw <- function(cross_sections, gauge_data, dx = 1, n.obs)
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
    if (dbf[seg] < 0){next} # error check
    
    # Make a vector of WSE(t) values to enter in with
    
    WSE_bf <- b.min[seg] + dbf[seg] # bankfull WSE
    WSE <- x1*WSE_bf

    # Estimate the flow width for each WSE value
    w <- get_width(WSE, x[[seg]], b[[seg]], delx = dx)
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
