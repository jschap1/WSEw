#' Calculate A0 for Pepsi Challenge data
#' 
#' Calculates A0, the flow area below the observed water surface elevation.
#' @param t time step
#' @param WSEw WSE-w data frame for a reach (not sorted)
#' @param A flow area for reach r at time t
#' @examples r <- 1
#' t <- 365
#' A0[r,t] <- calc_A0_pepsi(WSEw[[r]], t, A[r,t])

calc_A0_pepsi <- function(WSEw, t, A)
{
  
  hr <- WSEw$WSE
  wr <- WSEw$w
  WSEw.sorted <- data.frame(WSE = hr[order(hr)], w = wr[order(hr)])
  
  WSE1 <- WSEw.sorted$WSE[WSEw.sorted$w >= WSEw$w[t]]
  w1 <- WSEw.sorted$w[WSEw.sorted$w >= WSEw$w[t]]
  rWSEw <- data.frame(WSE = WSE1, w = w1)
  
  A.prime <- calc_overlying_A(rWSEw) # A = A0 + A.prime
  A0 = A - A.prime
  
  if (length(A0) == 0) # to account for when the highest observation is reached
  {
    A0 <- 0
  }
  
  return(A0)
}

# A0[1,77] # is negative, for example. Why?
# It may be because I did averaging. Perhaps I need to calculate A0 at the cross section scale, 
# then average to the reach scale. Try this out. (It doesn't help.)

# ------------------------------------------------------------------------------------------------

#' Calculates area overlying the lowest observed WSE
#' 
#' That is, calculates A.prime = A - A0
#' @param WSEw must be sorted, include only values down to the desired WSE level
#' @export
calc_overlying_A <- function(WSEw)
{
  
  w <- WSEw$w
  h <- WSEw$WSE
  
  A.prime <- 0
  nn <- length(w)
  for (i in 1:(nn-1))
  {
    A.prime <- A.prime + 0.5*(w[i]+w[i+1])*(h[i+1]-h[i])
  }
  
  return(A.prime)
}

