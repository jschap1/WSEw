#' Calculate WSE-w
#' 
#' Calculates WSE-w pairs for each cross section
#' @export
#' @param cross_sections list of x coordinates, bed elevations, and depths for each transect
#' @param interval number of points to sample in each cross-section. Unlike calc_WSE, interval does not determine the number of data points, which are instead constrained by the SWOT resolution (50m) and the river width
#' @param dx used for computing width with trapezoidal rule
#' @details Makes a vector of WSE values from the minimum bed elevation to the bankfull WSE
#' Computes flow width for each WSE-w value.
#' Puts the estimated flow widths into the WSEw data frame
#' @return WSEw data frame containing WSE and flow width values
#' @examples 
#' xWSEw <- calc_WSEw(cross_sections, interval = 0.05, dx = 1)

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
    if (dbf[seg] <= 0){next} # error check
    WSE <- seq(b.min[seg], b.min[seg] + dbf[seg], by = interval*dbf[seg])
    
    # Estimate the flow width for each WSE value
    w <- get_width(WSE, x[[seg]], b[[seg]], delx = dx)
    
    # Put WSE-w info into a data frame
    WSEw <- data.frame(WSE = WSE, w = w)

    xWSEw[[seg]] <- WSEw
  }
  
  return(xWSEw)
}

#' Remove null cross sections
#' 
#' Removes null cross sections, those with no measurements (optional)
#' @export
#' @param return_ind boolean telling the function whether or not to return rm.ind
rm_null_xs <- function(xWSEw, return_ind = FALSE)
{
  rm.ind <- which(unlist(lapply(xWSEw, is.null)))
  if (any(rm.ind))
  {
    xWSEw <- xWSEw[-rm.ind]
  }
  
  if (return_ind)
  {
    result <- list(WSEw = xWSEw, rm.ind = rm.ind)
    return(result)
  } else 
  {
    return(xWSEw)
  }
}


