#' Calculate WSEw from stage time series
#' 
#' Calculates water surface elevation-width pairs from stage time series to simulate the data available to SWOT
#' Corrupts the synthetic measurements, as well
#' @param cross_sections list of cross sections
#' @param h gauged water surface elevations
#' @param t vector of times corresponding to gauged water surface elevations
#' @param tobs vector of times at which to make an observation 
#' @param r cross section number
#' @examples 
#' @export
#' @details Similar to calc_WSEw, but selects WSE values to enter in with based on a time series of gauge measurements.
#' Also incorporates some of the functionality from observe().
#' This function is a work in progress.
 
calc_WSEw_from_stage_data <- function(cross_sections, h, t, tobs, r, sd_w = 0, sd_h = 0)
{
  
  x <- cross_sections$x[[r]]
  b <- cross_sections$b[[r]]
  d <- cross_sections$d[[r]]
  
  # Pre-processing
  dbf <- max(d, na.rm = TRUE) # bankfull depth (m)
  b.min <- min(b, na.rm = TRUE) # bottom depth
  
    
  if (dbf < 0) # error check
  {
    warning("dbf<0")
    return(NULL)
  } 

    
  #   WSE <- seq(b.min[seg], b.min[seg] + dbf[seg], by = interval*dbf[seg])
  #   
  #   # Estimate the flow width for each WSE value
  #   w <- get_width(WSE, x[[seg]], b[[seg]], delx = dx)
  #   
  #   # Put WSE-w info into a data frame
  #   WSEw <- data.frame(WSE = WSE, w = w)
  #   
  #   WSEw[[seg]] <- WSEw
  # 
  # 
  # # Remove null cross sections, those with no measurements (optional)
  # rm.ind <- which(unlist(lapply(WSEw, is.null)))
  # if (any(rm.ind))
  # {
  #   WSEw <- WSEw[-rm.ind] # this doesn't always seem to work...
  # }
  
  
  
  return(WSEw)
}