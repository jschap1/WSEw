#' Fit height-width data (parallel)
#'
#' Parallelized fitting functions for model fitting with foreach
#' @export
#' @param r cross section number
#' @details Saves fitted model (linear, slope break, multiple slope break, nonlinear, or nonlinear slope break) 
#' and synthetic height-width observations.
#' @return list of z0, A, W0, and A0 predictions
#' @example pred_vals <- pred_linear_par(r = 1)

# Make observations
observe_par <- function(r, WSEw)
{
  WSEw_obs <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels) # loop over exposure levels
  {
    WSEw_obs[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
      WSEw_obs[[k]][[m]] <- observe(WSEw = rWSEw[[r]], exposure = expo[k])
    }
  }
  # save data for each cross section (to avoid bulky files)
  obsname <- paste0("obs/WSEw_obs_r_", r, ".rds")
  saveRDS(WSEw_obs, file = file.path(exp_dir, obsname))
  return(0)
}

# -------------------------------------------------------------------------------------------------------------------------------------

#' @export
# Fit linear model
fit_linear_par <- function(r)
{
  obsname <- paste0("obs/WSEw_obs_r_", r, ".rds") # load observations for this reach
  WSEw_obs <- readRDS(file.path(exp_dir, obsname))
  lf <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels)
  {
    lf[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
        lf[[k]][[m]] <- fit_linear(WSEw_obs[[k]][[m]], h = 5)
    }
  }
  lf_name <- paste0("lf/lf_", "r_", r, ".rds")
  saveRDS(lf, file = file.path(exp_dir, lf_name))
  return(0)
}

# r <- 7
# lf_name <- paste0("lf/lf_", "r_", r, ".rds")
# lf <- readRDS(file.path(exp_dir, lf_name))
# lf[[1]][[10]]

#' @export
# Fit SB model
fit_sb_par <- function(r)
{
  obsname <- paste0("obs/WSEw_obs_r_", r, ".rds")
  WSEw_obs <- readRDS(file.path(exp_dir, obsname))
  sb <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels)
  {
    sb[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
        sb[[k]][[m]] <- fit_slopebreak(WSEw_obs[[k]][[m]], 
                                       multiple_breaks = FALSE, 
                                       continuity = TRUE, 
                                       minlen = 5)
    }
  }
  sb_name <- paste0("sb/sb_", "r_", r, ".rds")
  saveRDS(sb, file = file.path(exp_dir, sb_name))
  return(0)
}

#' @export
# Fit nonlinear model
fit_nl_par <- function(r)
{
  obsname <- paste0("obs/WSEw_obs_r_", r, ".rds") # load observations for this reach
  WSEw_obs <- readRDS(file.path(exp_dir, obsname))
  nl <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels)
  {
    nl[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
      model1 <- fit_nonlinear(WSEw_obs[[k]][[m]], h=5) # used to be in try() but i don't think it's nec anymore
      if (attributes(model1)$ef==0)
      {
        nl[[k]][[m]] <- list(z0 = as.numeric(coef(model1))[1],
             a = as.numeric(coef(model1))[2],
             s = as.numeric(coef(model1))[3],
             WSEbf = max(fitted(model1)),
             ef = attributes(model1)$ef)
      } else
      {
        nl[[k]][[m]] <- list(z0 = NA,
                             a = NA,
                             s = NA,
                             WSEbf = NA,
                             ef = attributes(model1)$ef)
      }
      # print(paste("finished processing for replicate", m, "of 500"))
    }
    print(paste("finished processing for exposure level", k, "of 19"))
  }
  
  nl_name <- paste0("nl/nl_", "r_", r, ".rds")
  saveRDS(nl, file = file.path(exp_dir, nl_name))
  return(0)
}

#' @export
# Fit NLSB model
fit_nlsb_par <- function(r)
{
  obsname <- paste0("obs/WSEw_obs_r_", r, ".rds") # load observations for this reach
  WSEw_obs <- readRDS(file.path(exp_dir, obsname))
  nlsb <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels)
  {
    nlsb[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
      model1 <- tryCatch( # note: it is imperative to have model1 <- tryCatch() to return output at all
        {
          model1 <- fit_nlsb(WSEw_obs[[k]][[m]], h = 10)
        }, 
        error = function(e) 
        {
          # print("error: nonlinear fit did not converge")
          model1 <- NULL
        }
      )
      if (attributes(model1)$ef==0)
      {
        nlsb[[k]][[m]] <- list(z0 = as.numeric(coef(model1[[1]]))[1],
                             a1 = as.numeric(coef(model1[[1]]))[2],
                             s1 = as.numeric(coef(model1[[1]]))[3],
                             a2 = as.numeric(coef(model1[[2]]))[1],
                             s2 = as.numeric(coef(model1[[2]]))[2],
                             WSEbf = max(fitted(model1[[2]])),
                             sb.ind = attributes(model1)$sb.ind,
                             ef = attributes(model1)$ef)
      } else
      {
        nlsb[[k]][[m]] <- list(z0 = NA,
                               a1 = NA,
                               s1 = NA,
                               a2 = NA,
                               s2 = NA,
                               WSEbf = NA,
                               sb.ind = NA,
                               ef = attributes(model1)$ef)
      }
    }
    print(paste("finished processing for exposure level", k, "of 19"))
  }
  nlsb_name <- paste0("nlsb/nlsb_", "r_", r, ".rds")
  saveRDS(nlsb, file = file.path(exp_dir, nlsb_name))
  return(0)
}