#' Plots bias vs. exposure level
#' Also calculates bias
#' 
#' @param expo vector of exposure levels
#' @param pred_vals predicted values of bottom elevation, flow area, wetted perimeter, etc. A list with 5 elements.
#' @param true_vals true values
#' @param na.rm
#' @param ... other arguments supported by plot
#' @export
#' @details 
#' Must be five pred_vals inputs, corresponding to the linear, slope break, multiple slope break, nonlinear, and nlsb methods
#' @return bias a data frame showing bias at each exposure level, averaged across all river reaches

plot_bias <- function(expo, pred_vals, true_vals, na.rm = TRUE, ...)
{
 
  n_levels <- dim(pred_vals[[1]])[2]
  
  l <- pred_vals[[1]]
  sb <- pred_vals[[2]]
  sbm <- pred_vals[[3]]
  nl <- pred_vals[[4]]
  nlsb <- pred_vals[[5]]
  
  # Storing in a data frame
  bias <- array(dim = c(n_levels, 5))
  for (k in 1:n_levels)
  {
    bias[k,1] <- mean(l[,k] - true_vals, na.rm = na.rm)
    bias[k,2] <- mean(sb[,k] - true_vals, na.rm = na.rm)
    bias[k,3] <- mean(sbm[,k] - true_vals, na.rm = na.rm)
    bias[k,4] <- mean(nl[,k] - true_vals, na.rm = na.rm)
    bias[k,5] <- mean(nlsb[,k] - true_vals, na.rm = na.rm)
  }
  bias <- as.data.frame(bias)
  names(bias) <- c("l","sb","sbm","nl","nlsb")
  
  # Plot bias vs. exposure level
  plot(100*expo, bias$l, col = "red", type = "l", xlab = "Channel exposure (%)", ...)
  lines(100*expo, bias$sb, col = "orange")
  lines(100*expo, bias$sbm, col = "purple")
  lines(100*expo, bias$nl, col = "green")
  lines(100*expo, bias$nlsb, col = "blue")
  abline(0,0)
  # legend("bottomright", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
  #        col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)
  
  return(bias)
   
}
