#' Computes RMSE
#' 
#' Computes RMSE between the model predictions of hydraulic parameters and their true values
#' @export
#' @param pred_vals predicted values of bottom elevation, flow area, wetted perimeter, etc. A list with 5 elements.
#' @param true_vals true values
#' @details 
#' Must be five pred_vals inputs, corresponding to the linear, slope break, multiple slope break, nonlinear, and nlsb methods
#' @return rmse, a data frame showing rmse at each exposure level, averaged across all river reaches
#' @example compute_rmse(z0.pred, z0.true.ra)

compute_rmse <- function(pred_vals, true_vals)
{
  
  n_levels <- dim(pred_vals[[1]])[2]
  
  l <- pred_vals[[1]]
  sb <- pred_vals[[2]]
  sbm <- pred_vals[[3]]
  nl <- pred_vals[[4]]
  nlsb <- pred_vals[[5]]
  
  # Storing in a data frame
  rmse <- array(dim = c(n_levels, 5))
  for (k in 1:n_levels)
  {
    rmse[k,1] <- ((1/length(true_vals))*t(l[,k] - true_vals)%*%(l[,k] - true_vals))^0.5
    rmse[k,2] <- ((1/length(true_vals))*t(sb[,k] - true_vals)%*%(sb[,k] - true_vals))^0.5
    rmse[k,3] <- ((1/length(true_vals))*t(sbm[,k] - true_vals)%*%(sbm[,k] - true_vals))^0.5
    rmse[k,4] <- ((1/length(true_vals))*t(nl[,k] - true_vals)%*%(nl[,k] - true_vals))^0.5
    rmse[k,5] <- ((1/length(true_vals))*t(nlsb[,k] - true_vals)%*%(nlsb[,k] - true_vals))^0.5
  }
  rmse <- as.data.frame(rmse)
  names(rmse) <- c("l","sb","sbm","nl","nlsb")
  
  return(rmse)
  
}

# mse <- (1/length(true_vals))*t(l[,k] - true_vals)%*%(l[,k] - true_vals)
# rmse <- mse^0.5