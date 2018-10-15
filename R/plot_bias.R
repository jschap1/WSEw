#' Plots bias vs. exposure level
#' 
#' Also calculates bias
#' @export
#' @param expo vector of exposure levels
#' @param pred_vals predicted values of bottom elevation, flow area, wetted perimeter, etc. A list with 5 elements.
#' @param true_vals true values
#' @param na.rm
#' @param ... other arguments supported by plot
#' @details 
#' Must be five pred_vals inputs, corresponding to the linear, slope break, multiple slope break, nonlinear, and nlsb methods
#' @return bias a data frame showing bias at each exposure level, averaged across all river reaches
#' @examples
#' pred_z0 <- list(z0.l, z0.sb, z0.nl, z0.nlsb)
#' plot_bias(expo, pred_z0, z0.true, na.rm = TRUE, 
#'          main = "z0 bias vs. exposure level, no meas. error", ylab = "Bias (m)", ylim = c(-2,2))

plot_bias <- function(expo, pred_vals, true_vals, na.rm = TRUE, ...)
{
 
  n_levels <- dim(pred_vals[[1]])[2]
  
  l <- pred_vals[[1]]
  sb <- pred_vals[[2]]
  nl <- pred_vals[[3]]
  nlsb <- pred_vals[[4]]
  
  # # Handle NaN values
  # l[is.na(l)] <- 9999
  # sb[is.na(sb)] <- 9999
  # nl[is.na(nl)] <- 9999
  # nlsb[is.na(nlsb)] <- 9999
  # true_vals[is.na(true_vals)] <- 9999
  
  # Storing in a data frame
  bias <- array(dim = c(n_levels, 4))
  for (k in 1:n_levels)
  {
    bias[k,1] <- mean(l[,k] - true_vals, na.rm = na.rm)
    bias[k,2] <- mean(sb[,k] - true_vals, na.rm = na.rm)
    bias[k,3] <- mean(nl[,k] - true_vals, na.rm = na.rm)
    bias[k,4] <- mean(nlsb[,k] - true_vals, na.rm = na.rm)
  }
  bias <- as.data.frame(bias)
  names(bias) <- c("l","sb","nl","nlsb")
  
  # Plot bias vs. exposure level
  plot(100*expo, bias$l, col = "red", type = "l", xlab = "Channel exposure (%)", ...)
  lines(100*expo, bias$sb, col = "orange")
  lines(100*expo, bias$nl, col = "green")
  lines(100*expo, bias$nlsb, col = "blue")
  abline(0,0)
  
  return(bias)
   
}
