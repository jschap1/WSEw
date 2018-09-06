#' Calculate A0 variance
#' 
#' Calculates A0 variance for the linear case
#' @export
#' @param m linear model
#' @param sd_h standard deviation of WSE measurement error (m)
#' @details 
#' Assumes no width measurement error
#' Assumes that WSE measurement errors are independent of WSE
#' 
#' @example calc_A0_variance(lf[[1]][[17]])

calc_A0_variance <- function(m, sd_h = 0.1)
{
  w1 <- min(m$model$w) # minimum observed width
  var_z0 <- vcov(m)[1,1] # variance of the parameter estimate
  varA0 <- (1/4)*w1^2*(sd_h^2 + var_z0)
}

