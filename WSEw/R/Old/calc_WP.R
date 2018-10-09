#' Calculate Wetted Perimeter
#' 
#' Calculates wetted perimeter for a particular water surface elevation, ranging from empty channel to bankfull conditions
#' @param x distance from horizontal datum
#' @param b channel bed elevation
#' @param WSE water surface elevation
#' @param method "neal" for Jeff Neal's approximation in Neal et al. (2015), "wide" for wide-channel approximation
#' @importFrom pracma trapz
#' @details
#' Uses a numerical method to approximate the integral of the difference between WSE and bed elevation
#' The numerical method only counts area where the WSE is higher than the bed elevation
#' Needs more work to figure out how to implement arbitrary spatial discretization
#' Performs poorly for a small number of data points or WSE close to bed elevation, see the example
#' @examples calc_WP(x = cross_sections$x[[1]], b = cross_sections$b[[1]], WSE = 138)
#' x <- c(1,1,1,2,3,4,4,4)
#' b <- c(3,2,1,1,1,1,2,3)
#' WP <- calc_WP(x, b, WSE = 3)
#' WP <- calc_WP(x, b, WSE = 2)

calc_WP <- function(x, b, WSE, method = "wide")
{
  
  if (method == "wide")
  {
    w <- get_width(WSE, x, b, delx = 0.1) # find the width, fails if WSE > bankfull
    d <- WSE- min(b) # find the depth
    # evaluate whether the wide-channel assumption is appropriate
    if (w==0)
    {
      WP <- 0
      return(WP)
    }
    if (d/w > 0.05) # check literature (Dingman 2009) for a good threshold
    {
      warning("d/w > 0.05. Using wide-channel assumption is likely to introduce significant error for this channel.")
    }
    WP <- w
    
  }else if (method == "neal")
  {
    print("full method not implemented yet")
    WP <- NA
    
    
  }else if (method == "full")
  {
    print("full method not implemented yet")
    WP <- NA
  }
  
  return(WP)
  
}
