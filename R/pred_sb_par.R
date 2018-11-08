#' Slope break
#' 
#' Applies fit_slopebreak in parallel across multiple processors on one machine
#' @export
#' @param r cross section number
#' @param exclude boolean that determines whether or not to exclude predictions from physically unrealistic model fits
#' @details Loads fitted model (linear, slope break, multiple slope break, nonlinear, or nonlinear slope break)
#' and uses it to predict hydraulic parameters
#' error.flag values: 0 = no error, 1 = negative slope, 3 = no model was fit
#' @return list of z0, A, W0, and A0 predictions

pred_sb_par <- function(r, exclude = FALSE)
{
  sb_name <- paste0("sb/sb_", "r_", r, ".rds")
  sb <- readRDS(file.path(exp_dir, sb_name))
  z0.sb <- array(dim = c(n_exp_levels, M))
  A.sb <- array(dim = c(n_exp_levels, M))
  WP.sb <- array(dim = c(n_exp_levels, M))
  A0.sb <- array(dim = c(n_exp_levels, M))
  error.flag <- array(dim = c(n_exp_levels, M), data = 0)
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      if (attributes(sb[[k]][[m]])$ef==0) 
      {
        
        # Check that physical-realism constraints are satisfied
        # Must have non-negative slope
        a1 <- coef(sb[[k]][[m]][[1]])[2] # slope coefficient
        if (a1 < 0)
        {
          error.flag[k,m] <- 1
          if (exclude) # skip bc everything is noise
          {
            next
          } else # try to estimate parameters anyway
          {
            m1 <- sb[[k]][[m]][[1]] # lower model
            m2 <- sb[[k]][[m]][[2]] # upper model 
            wsb <- max(m1$model$w)
            wbf <- max(m2$model$`I(w - wb)`) + wsb
            WSEsb <- max(m1$model$WSE) # choose max observed WSE
            WSEbf <- max(m2$model$WSE)
            z0.sb[k,m] <- min(sb[[k]][[m]][[1]]$model$WSE) # choose min observed WSE
            A.sb[k,m] <- sb_area(wsb, WSEsb, z0.sb[k,m], WSEbf, wbf) 
            WP.sb[k,m] <- calc_WP_sb(wbf, WSEbf, z0.sb[k,m], wsb, WSEsb)
            A0.sb[k,m] <- 0
          }
        } else
        {
          z0.sb[k,m] <- as.numeric(coef(sb[[k]][[m]][[1]])[1])
          A.sb[k,m] <- calc_model_A(sb[[k]][[m]], type = "sb")
          WP.sb[k,m] <- calc_model_WP(sb[[k]][[m]], type = "sb")
          A0.sb[k,m] <- calc_model_A0(sb[[k]][[m]], type = "sb", pos.only = FALSE)
        }
      }
      else # no model was fit
      {
        return(list(z0 = NA, A = NA, WP = NA, A0 = NA, ef = 3))
      }
    }
  }
  return(list(z0 = z0.sb, A = A.sb, WP = WP.sb, A0 = A0.sb, ef = error.flag))
}