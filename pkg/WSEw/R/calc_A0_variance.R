#' Calculate variance of the A0 estimate
#' 
#' Currently only handles linear models, and also it assumes there is no width measurement error

calc_A0_variance <- function(model, type = "linear", sd_wse = 0.1)
{
  if (type == "linear")
  {
    var.z0.est <- vcov(model)[1,1] # variance of z0 parameter estimate
    w1 <- min(model$model$w) # lowest observed width
    v1 <- sd_wse^2 + var.z0.est # variance of z0 prediction error
    v2 <- (w1^2/4)*(v1) # variance of A0 prediction error
    sd.A0.l <- v2^0.5 # standard deviation of A0 prediction error
  }
  return(v2) # returns the variance
}
