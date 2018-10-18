#' Thin WSEw
#' 
#' Thins out the WSEw data if they are too dense
#' @param WSEw WSEw
#' @param n number of data points in each of the new WSEw[[r]] data frames
#' @example tWSEw <- thin_WSEw(WSEw, n)
thinWSEw <- function(WSEw, n)
{
  nr <- length(WSEw)
  tWSEw <- vector(length = nr, "list")
  for (r in 1:nr)
  {
    tWSEw[[r]] <- data.frame(WSE = approx(WSEw[[r]]$WSE, n = n)$y,
                             w = approx(WSEw[[r]]$w, n = n)$y)
  }
  return(tWSEw)
}