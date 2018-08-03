#' Get maxdiff
#' 
#' Gets the maximum difference in the slope between any two points in a series
#' @param WSE WSE
#' @param w width
#' @export 

get_maxdiff <- function(WSE, w)
{
  # For use with Mersel's linear method
  
  nn <- length(WSE)
  s1 <- vector(length = nn)
  for (k in 1:(nn-1))
  {
    s1[k] <- (WSE[k+1]-WSE[k])/(w[k+1]-w[k]) # forward difference
  }
  s1[nn] <- s1[nn-1]
  # remove infinite values (where width did not change)
  s1[s1==Inf] <- NA
  s1[s1==-Inf] <- NA
  
  maxdiff <- max(s1, na.rm = TRUE) - min(s1, na.rm=TRUE)
  
  return(maxdiff)
}