#' Thin WSEw
#' 
#' Thins out the WSEw data if they are too dense
#' @param WSEw WSEw
#' @param n number of data points in each of the new WSEw[[r]] data frames
#' @example tWSEw <- thinWSEw(WSEw, n)
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

# ------------------------------------------------------------------------------------------------
#' Sample WSEw
#' 
#' Samples the WSEw data every n timesteps
#' @param WSEw WSEw
#' @param n number of time steps between successive observations
#' @param first_overpass time step of first overpass
#' @example tWSEw <- sampleWSEw(WSEw, n)
sampleWSEw <- function(WSEw, n, first_overpass = 1)
{
  nr <- length(WSEw)
  nt <- length(WSEw[[1]]$w)
  tWSEw <- vector(length = nr, "list")
  for (r in 1:nr)
  {
    ind <- seq(first_overpass, nt, by = n)
    tWSEw[[r]] <- data.frame(WSE = WSEw[[r]]$WSE[ind], w = WSEw[[r]]$w[ind])
  }
  return(tWSEw)
}

