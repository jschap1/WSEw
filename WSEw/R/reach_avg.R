#' Reach Average
#'
#' Performs reach averaging on the WSE, w data
#' @export
#' @param l = reach length, m. Default is 10 km
#' @param res = cross-section spacing, m
#' @param xWSEw WSEw pairs for each cross section: xWSEw
#' @details Tailor made for WSEw data returned by the calc_WSEw function.
#' @return rWSEw Reach averaged WSE and w pairs
#' @examples rWSEw <- reach_avg(xWSEw, l = 10000, res = 50)

reach_avg <- function(xWSEw, l = 10000, res = 5)
{

  pix <- l/res # the number of sections in a reach
  nseg <- length(xWSEw) # number of cross sections
  nn <- dim(xWSEw[[1]])[1]
  
  # Create matrices with WSE, w values
  WSEmat <- array(dim = c(nn, nseg))
  wmat <- array(dim = c(nn, nseg))
  for (seg in 1:nseg)
  {
    WSEmat[,seg] <- xWSEw[[seg]]$WSE
    wmat[,seg] <- xWSEw[[seg]]$w
  }
  
  # Do reach averaging
  nr <- nseg - (pix - 1) # number of reaches
  rWSEw <- vector(length = nr, "list")
  for (r in 1:nr)
  {
    rWSEw[[r]] <- data.frame(WSE = rowSums(WSEmat[,r:(r+(pix-1))], na.rm = TRUE)/pix, 
                             w = rowSums(wmat[,r:(r+(pix-1))], na.rm = TRUE)/pix)
  }
  
  # Developing the averaging rule
  # rWSEw[[r]] <- data.frame(WSE = rowSums(WSEmat[,1:4])/pix # may want to handle NAs here
  #                         , w = rowSums(wmat[,1:4])/pix)
  
  # Note: the wbf values as used for defining channel exposure have changed. 
  # Now, it would make sense to use max(reach average w) as the wbf value.
  
  return(rWSEw)
}

# --------------------------------------------------------------------------------------------------

#' Reach Average General
#' 
#' Reach average any quantity
#' @export
#' @param x quantity, a vector or a list of numbers
#' @param n number of segments making up the reach
#' @details If the reach-average length is 10 km and each individual segment is 5 m, then n = 2000
#' @return r.out reach averaged quantity
#' @examples
#' x <- c(1,2,3,4,4,6)
#' n <- 3
#' r <- ra(x,n)

ra <- function(x,n)
{
  nseg <- length(x) # number of cross sections
  nr <- nseg - (n - 1) # number of reaches
  r.out <- vector(length = nr)
  for (r in 1:nr)
  {
    r.out[r] <- mean(x[r:(r+(n-1))])
  }
  
  return(r.out)
}

# --------------------------------------------------------------------------------------------------


