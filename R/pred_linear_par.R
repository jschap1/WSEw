#' Predict parameters (parallel)
#'
#' Parallelized functions for main_d6 hydraulic parameter predictions in parallel with foreach
#' for WSE-w fits.
#' @export
#' @param r cross section number
#' @param exclude boolean that determines whether or not to exclude predictions from physically unrealistic model fits
#' @details Loads fitted model (linear, slope break, multiple slope break, nonlinear, or nonlinear slope break)
#' and uses it to predict hydraulic parameters
#' error.flag values: 0 = no error, 1 = negative slope
#' @return list of z0, A, W0, and A0 predictions
#' @example pred_vals <- pred_linear_par(r = 1)

# Linear
pred_linear_par <- function(r, exclude = FALSE)
{
  z0.l <- array(dim = c(n_exp_levels, M)) # initialize
  A.l <- array(dim = c(n_exp_levels, M))
  WP.l <- array(dim = c(n_exp_levels, M))
  A0.l <- array(dim = c(n_exp_levels, M))
  error.flag <- array(dim = c(n_exp_levels, M), data = 0)
  lf_name <- paste0("lf/lf_", "r_", r, ".rds") # load fitted model
  lf <- readRDS(file.path(exp_dir, lf_name))
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      # if statements handle cases where the model is NULL/no model was fit
      if (class(lf[[k]][[m]]) == "lm")
      {
        
        # Check that physical-realism constraints are satisfied
        # Must have non-negative slope
        a <- coef(lf[[k]][[m]])[2] # slope coefficient
        if (a < 0)
        {
          error.flag[k,m] <- 1
          if (exclude) # skip bc everything is noise
          {
            next
          } else # attempt to estimate parameters anyway
          {
            wbf <- max(lf[[k]][[m]]$model$w)
            WSEbf <- max(lf[[k]][[m]]$model$WSE) # choose max observed WSE
            z0.l[k,m] <- min(lf[[k]][[m]]$model$WSE) # choose min observed WSE
            A.l[k,m] <- 0.5*wbf*(WSEbf - z0.l[k,m])
            WP.l[k,m] <- calc_WP_linear(wbf, WSEbf, z0.l[k,m])
            A0.l[k,m] <- 0 # since we've assumed z0 = min(h_obs)
          }
        } else
        {
          z0.l[k,m] <- as.numeric(coef(lf[[k]][[m]])[1])
          A.l[k,m] <- calc_model_A(lf[[k]][[m]], type = "linear")
          WP.l[k,m] <- calc_model_WP(lf[[k]][[m]], type = "linear")
          A0.l[k,m] <- calc_model_A0(lf[[k]][[m]], type = "linear", pos.only = FALSE)
        }
        
      }
    }
  }
  return(list(z0 = z0.l, A = A.l, WP = WP.l, A0 = A0.l, ef = error.flag))
}
