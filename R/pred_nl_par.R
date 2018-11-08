#' Nonlinear
#' 
#' Applies fit_nonlinear in parallel across multiple processors on one machine
#' @export
#' @param r cross section number
#' @param exclude boolean that determines whether or not to exclude predictions from physically unrealistic model fits
#' @param w1 minimum observed width, nr by n_exp_level by M array
#' @param h1 minimum observed height, nr by n_exp_level by M array
#' @details Loads fitted model (linear, slope break, multiple slope break, nonlinear, or nonlinear slope break)
#' and uses it to predict hydraulic parameters
#' error.flag values: 0 = no error, 1 = negative slope, 2 = negative concavity, 3 = no model was fit
#' There is a problem when s is negative because the predictions at w = 0 are -Inf
#' When there is negative concavity, the linear model is used to make predictions, instead of the nonlinear model
#' @return list of z0, A, W0, and A0 predictions

pred_nl_par <- function(r, WSEw, w1, h1, exclude = FALSE)
{
  # Load fitted NL model
  nl_name <- paste0("nl/nl_", "r_", r, ".rds")
  if (file.exists(file.path(exp_dir, nl_name))) # error check in case no model was fit for this cross section
  {
    nl <- readRDS(file.path(exp_dir, nl_name))
  } else
  {
    return(list(z0 = NA, A = NA, WP = NA, A0 = NA, ef = 3))
  }
  
  # Load fitted L model
  lf_name <- paste0("lf/lf_", "r_", r, ".rds")
  lf <- readRDS(file.path(exp_dir, lf_name))
  
  # ----------------------------------------------------------------------------------------------------
  
  # Initialize
  z0.nl <- array(dim = c(n_exp_levels, M))
  A.nl <- array(dim = c(n_exp_levels, M))
  WP.nl <- array(dim = c(n_exp_levels, M))
  A0.nl <- array(dim = c(n_exp_levels, M))
  error.flag <- array(dim = c(n_exp_levels, M), data = 0)
  
  # Loop over exposure levels and Monte Carlo replicates
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      if (all(!is.na(nl[[k]][[m]]))) # if an NL model has been fit
      {
        
        # Check that physical-realism constraints are satisfied
        # Must have non-negative slope and cannot be concave down

        a <- nl[[k]][[m]]$a # slope coefficient
        s <- nl[[k]][[m]]$s # shape parameter
        
        if (a>=0 & s>=1) # physically realistic case
        {
          # positive slope, concave up
          z0.nl[k,m] <- nl[[k]][[m]]$z0
          A.nl[k,m] <- calc_model_A(nl[[k]][[m]], type = "nl", WSEw = WSEw[[r]])
          WP.nl[k,m] <- calc_model_WP(nl[[k]][[m]], type = "nl", w = WSEw[[r]]$w)
          A0.nl[k,m] <- calc_model_A0(nl[[k]][[m]], type = "nl", w1 = w1[r,k,m], h1 = h1[r,k,m], pos.only = FALSE)
        } 
        
        else if ((a<=0 & s>=1) | (a<=0 & s<=1)) # negative slope case
        {
          error.flag[k,m] <- 1
          if (exclude) # skip bc everything is noise
          {
            next
          } else # try to estimate parameters anyway
          {
            # negative slope, assume z0 = min(h_obs), linear h-w relationship
            warning("exclude == FALSE functionality has not been added yet for nonlinear models")
            z0.nl[k,m] <- NA
            A.nl[k,m] <- NA
            WP.nl[k,m] <- NA
            A0.nl[k,m] <- NA
          }
        } 
        
        else if (a>=0 & s<=1) # positive slope, but negative concavity case
        {
          error.flag[k,m] <- 2
          warning("Nonlinear fit is concave down, switching to linear fit")
          if (class(lf[[k]][[m]]) == "lm")
          {
            a <- coef(lf[[k]][[m]])[2]
            if (a < 0) # if linear model has negative slope
            {
              if (exclude) # if excluding nonrealistic cases, returns NA
              {
                next
              } else # otherwise, predicts parameters as best it can (not justified, don't do it)
              {
                wbf <- max(lf[[k]][[m]]$model$w)
                WSEbf <- max(lf[[k]][[m]]$model$WSE) # choose max observed WSE
                z0.nl[k,m] <- min(lf[[k]][[m]]$model$WSE) # choose min observed WSE
                A.nl[k,m] <- 0.5*wbf*(WSEbf - z0.l[k,m])
                WP.nl[k,m] <- calc_WP_linear(wbf, WSEbf, z0.l[k,m])
                A0.nl[k,m] <- 0 # since we've assumed z0 = min(h_obs)
              }
            } 
            else # if linear model has positive slope
            {
              z0.nl[k,m] <- as.numeric(coef(lf[[k]][[m]])[1])
              A.nl[k,m] <- calc_model_A(lf[[k]][[m]], type = "linear")
              WP.nl[k,m] <- calc_model_WP(lf[[k]][[m]], type = "linear")
              A0.nl[k,m] <- calc_model_A0(lf[[k]][[m]], type = "linear", pos.only = FALSE)
            }
          }
        }
        
        
      }
    }
  }
  return(list(z0 = z0.nl, A = A.nl, WP = WP.nl, A0 = A0.nl, ef = error.flag))
}