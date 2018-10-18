#' Make observations
#' 
#' Sets the exposure level and adds error to the observations.
#' @param WSEw
#' @param exposure exposure level as a fraction of bankfull depth
#' @param sd_wse standard deviation of expected WSE error
#' @param sd_w standard deviation of expected width error
#' @export

observe <- function(WSEw, exposure, sd_wse = 0.1, sd_w = 10)
{
  
  # Keep observations based on the amount of exposed riverbank
  wbf <- max(WSEw$w) # bankfull width
  observed.ind <- which(WSEw$w/wbf > (1-exposure))
  WSEw_obs <- WSEw[observed.ind,]
  
  # corrupt observations
  nn <- dim(WSEw_obs)[1] # number of observations
  e_WSE <- rnorm(nn, mean = 0, sd = sd_wse)
  e_w <- rnorm(nn, mean = 0, sd = sd_w)
  
  WSEw_corr <- WSEw_obs
  WSEw_corr$WSE <- WSEw_obs$WSE + e_WSE
  WSEw_corr$w <- WSEw_obs$w + e_w
  
  # Constrain observations to be positive:
  neg.ind <- which(WSEw_corr<0, arr.ind = TRUE) 
  WSEw_corr[neg.ind] <- 0

  # Remove observations where width is less than 100 m
  narrow.ind <- which(WSEw_corr<=100, arr.ind = TRUE)[1]
  if (!is.na(narrow.ind))
  {
    WSEw_corr <- WSEw_corr[-narrow.ind,]
  }
  return(WSEw_corr)
}

