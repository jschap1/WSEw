#' Get Width
#' 
#' Numerically estimates the top width
#' @param WSE WSE values with which to enter into the WSE-w function
#' @param f1 function approximating the WSE-w relationship (like a spline)
#' @param x.max 
#' @param delx Default is 1.
#' @export
#' @examples w <- get_width(WSE, f1, x.max = max(cross.section$x), delx = dx)

get_width <- function(WSE, f1, x.max, delx = 1)
{
  
  n <- length(WSE)
  flow_width <- vector(length = n)
  
  b.trap <- f1(seq(0, x.max, by = delx))
  
  for (k in 1:n)
  {
    flow_width[k] <- delx*sum(WSE[k]>b.trap)
    
    # Once the WSE is higher than any other part of the bed elevation,
    # the width stops changing, as rectangular walls as assumed
    if (sum(WSE[k] > max(b.trap)))
    {
      flow_width[k] <- flow_width[k-1]
    }
  }
  
  return(flow_width)
}