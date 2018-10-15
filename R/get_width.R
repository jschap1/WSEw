#' Get Width
#' 
#' Numerically estimates the top width
#' @param WSE WSE values with which to enter into the WSE-w function
#' @param f1 function approximating the WSE-w relationship (like a spline)
#' @param x.max 
#' @param delx Default is 1. Accuracy improves with smaller delx
#' @export
#' @examples w <- get_width(WSE, f1, x.max = max(cross.section$x), delx = dx)
#' @details 

get_width <- function(WSE, x, b, delx = 1)
{
  
  x.max = max(x)
  n <- length(WSE)
  flow_width <- vector(length = n)
  
  # Fit a linear spline to the observed channel geometry
  f1 <- approxfun(x = x, y = b, method="linear", 
                  yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)
  b.est <- f1(seq(0, x.max, by = delx)) # compute interpolated b values
  
  for (k in 1:n)
  {
    flow_width[k] <- delx*sum(WSE[k]>b.est)
    
    # Once the WSE is higher than any other part of the bed elevation,
    # the width stops changing, as rectangular walls as assumed
    if (sum(WSE[k] > max(b.est)) & n>1)
    {
      flow_width[k] <- flow_width[k-1]
    }
  }
  
  return(flow_width)
}
