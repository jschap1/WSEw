#' Nonlinear slope break
#' 
#' Applies fit_nlsb in parallel across multiple processors on one machine
#' @export
#' @param r cross section number
#' @param exclude boolean that determines whether or not to exclude predictions from physically unrealistic model fits
#' @param w1 minimum observed width, nr by n_exp_level by M array
#' @param h1 minimum observed height, nr by n_exp_level by M array
#' @details Loads fitted model (linear, slope break, multiple slope break, nonlinear, or nonlinear slope break)
#' and uses it to predict hydraulic parameters
#' @return list of z0, A, W0, and A0 predictions

pred_nlsb_par <- function(r, WSEw, w1, h1, exclude = FALSE)
{
  
  # Load fitted NLSB model
  nlsb_name <- paste0("nlsb/nlsb_", "r_", r, ".rds")
  if (file.exists(file.path(exp_dir, nlsb_name))) # error check in case no model was fit for this cross section
  {
    nlsb <- readRDS(file.path(exp_dir, nlsb_name))
  } else
  {
    return(list(z0 = NA, A = NA, WP = NA, A0 = NA))
  }
  
  # Load fitted SB model
  sb_name <- paste0("sb/sb_", "r_", r, ".rds")
  sb <- readRDS(file.path(exp_dir, sb_name))
  
  # Initialize
  z0.nlsb <- array(dim = c(n_exp_levels, M))
  A.nlsb <- array(dim = c(n_exp_levels, M))
  WP.nlsb <- array(dim = c(n_exp_levels, M))
  A0.nlsb <- array(dim = c(n_exp_levels, M))
  error.flag <- array(dim = c(n_exp_levels, M), data = 0)
  
  # Loop over exposure levels and Monte Carlo replicates
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      if (all(!is.na(nlsb[[k]][[m]])))
      {
        
        # Check for physically realistic fits
        a1 <- nlsb[[k]][[m]]$a1 # slope coefficient
        s1 <- nlsb[[k]][[m]]$s1 # shape parameter
        
        if (a1>=0 & s1>=1) # physically realistic case
        {
          z0.nlsb[k,m] <- nlsb[[k]][[m]]$z0
          A.nlsb[k,m] <- calc_model_A(nlsb[[k]][[m]], type = "nlsb", WSEw = WSEw[[r]]) # there may be a bug in the type = nlsb code here
          WP.nlsb[k,m] <- calc_model_WP(nlsb[[k]][[m]], type = "nlsb", w = WSEw[[r]]$w)
          A0.nlsb[k,m] <- calc_model_A0(nlsb[[k]][[m]], type = "nlsb", w1 = w1[r,k,m], h1 = h1[r,k,m], pos.only = FALSE)
        } 
        
        else if ((a1<=0 & s1>=1) | (a1<=0 & s1<=1)) # negative slope case
        {
          error.flag[k,m] <- 1
          if (exclude) # skip bc everything is noise
          {
            next
          } else # try to estimate parameters anyway
          {
            # negative slope, assume z0 = min(h_obs), linear h-w relationship
            warning("exclude == FALSE functionality has not been added yet for nonlinear models")
            next
          }
        } 
        
        else if (a1>=0 & s1<=1) # positive slope, but negative concavity case
        {
          # positive slope, concave down, use SB fit
          error.flag[k,m] <- 2
          warning("Nonlinear fit is concave down, switching to slope break")
          if (class(sb[[k]][[m]][[1]]) == "lm")
          {
            a1 <- coef(sb[[k]][[m]][[1]])[2] # slope coefficient
            if (a1 < 0)
            {
              if (exclude) # skip bc everything is noise
              {
                next
              } else # try to estimate parameters anyway (no justification for this, don't do it)
              {
                m1 <- sb[[k]][[m]][[1]] # lower model
                m2 <- sb[[k]][[m]][[2]] # upper model 
                wsb <- max(m1$model$w)
                wbf <- max(m2$model$`I(w - wb)`) + wsb
                WSEsb <- max(m1$model$WSE) # choose max observed WSE
                WSEbf <- max(m2$model$WSE)
                z0.nlsb[k,m] <- min(sb[[k]][[m]][[1]]$model$WSE) # choose min observed WSE
                A.nlsb[k,m] <- sb_area(wsb, WSEsb, z0.sb[k,m], WSEbf, wbf) 
                WP.nlsb[k,m] <- calc_WP_sb(wbf, WSEbf, z0.sb[k,m], wsb, WSEsb)
                A0.nlsb[k,m] <- 0
              }
            } else # if SB model has positive slope
            {
              z0.nlsb[k,m] <- as.numeric(coef(sb[[k]][[m]][[1]])[1])
              A.nlsb[k,m] <- calc_model_A(sb[[k]][[m]], type = "sb")
              WP.nlsb[k,m] <- calc_model_WP(sb[[k]][[m]], type = "sb")
              A0.nlsb[k,m] <- calc_model_A0(sb[[k]][[m]], type = "sb", pos.only = FALSE)
            }
          }
        }
        
      }
    }
  }
  return(list(z0 = z0.nlsb, A = A.nlsb, WP = WP.nlsb, A0 = A0.nlsb, ef = error.flag))
}