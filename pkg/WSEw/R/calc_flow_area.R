#' Calculate Flow Area
#' 
#' Calculates cross sectional area of flow for a particular water surface elevation, ranging from empty channel to bankfull conditions
#' @param x distance from horizontal datum
#' @param b channel bed elevation
#' @param N spatial discretization, do not modify, except possibly for Riemann method
#' @param WSE water surface elevation
#' @param method "trap" "riemann" or "simpson" Default is trapezoidal.
#' @importFrom pracma trapz
#' @details
#' Uses a numerical method to approximate the integral of the difference between WSE and bed elevation
#' The numerical method only counts area where the WSE is higher than the bed elevation
#' Needs more work to figure out how to implement arbitrary spatial discretization
#' Need to do code testing, too. Confirmed that riemann and trap methods give similar results, but are they correct?
#' @example calc_A(x = cross_sections$x[[1]], b = cross_sections$b[[1]], WSE = 138)

calc_A <- function(x, b, WSE, N = length(x), method = "trap")
{
  
  if (method == "trap")
  {
    
    nn <- length(b)
    y <- rep(WSE, nn) - b # get difference btw WSE and bed elevation in locations where WSE>b(x)
    y[which(y<0)] <- 0
    A <- trapz(x, y)
    
  } else if (method == "riemann")
  {
    # Riemann sum (trapezoidal rule would be preferable)
      N <- length(x)
      delx <- (max(x)-min(x))/N
      A <- 0
      for (i in 1:N)
      {
        A <- A + max((WSE-b[i])*delx, 0) 
      }
  } else if (method == "simpson")
  {
    A <- NA
    print("simpson's method not yet implemented")
  }
  
  return(A)
}
