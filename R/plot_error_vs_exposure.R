#' Plot error vs. exposure level
#' 
#' Plots error vs. exposure level for a given method (L, SB, NL, NLSB)
#' @param err z0 error array
#' @param ... additional arguments to plot()
#' @example plot_error(z0.error.l, main = "title", ylab = "error (m)")
#' @export
#' @return mean.error, the mean vector
plot_error_vs_exposure <- function(err, ...)
{
  nr <- dim(err)[1]
  n_exp_levels <- dim(err)[2]
  mean.error <- array(dim = c(nr, n_exp_levels))
  for (r in 1:nr)
  {
    for (k in 1:n_exp_levels)
    {
      mean.error[r,k] <- mean(err[r,k,], na.rm = TRUE)
    }
  }
  
  matplot(100*seq(0.05, 0.95, by = 0.05), 
          t(mean.error), 
          type = "l", 
          lty = 1,
          col = rainbow(11),
          xlab = "Channel exposure (%)",
          ...
  ) 
  return(mean.error)
}

#' Plot breakpoint locations
#'
#' Plots breakpoint locations
#' @param brk_locs breakpoint locations
#' @param expo vector of exposure levels
#' @param mean.error mean error vector returned by plot_error_vs_exposure()
#' @param ... additional arguments for plot
#' @export
plot_brkpt_locs <- function(brk_locs, expo, mean.error, ...)
{
  nr <- length(brk_locs)
  for (r in 1:nr)
  {
    expo.ind <- which.min(abs(expo - brk_locs[r]))
    points(100*expo[expo.ind], mean.error[r,expo.ind], ...)
  }
}

