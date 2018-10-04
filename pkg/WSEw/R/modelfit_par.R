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
observe_par <- function(r)
{
  WSEw_obs <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels) # loop over exposure levels
  {
    WSEw_obs[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
      WSEw_obs[[k]][[m]] <- observe(WSEw = rWSEw[[r]], exposure = expo[k],
                                    sd_wse = 1e-6, sd_w = 1e-6)
      # WSEw_obs[[k]][[m]] <- observe(WSEw = xWSEw[[r]], exposure = expo[k])
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
  # Only run the code when there are at least 2 observations
  nn <- array(dim = c(n_exp_levels, M))
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      nn[k,m] <- dim(WSEw_obs[[k]][[m]])[1]
    }
  }
  lf <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels)
  {
    lf[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
      if (nn[k,m]>2)
      {
        lf[[k]][[m]] <- fit_linear(WSEw_obs[[k]][[m]])
      }
    }
  }
  lf_name <- paste0("lf/lf_", "r_", r, "_test.rds")
  saveRDS(lf, file = file.path(exp_dir, lf_name))
  return(0)
}

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
      try(sb[[k]][[m]] <- fit_slopebreak(WSEw_obs[[k]][[m]], multiple_breaks = FALSE, continuity = TRUE, minlen = 3))
    }
  }
  sb_name <- paste0("sb/sb_", "r_", r, "_test.rds")
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
      try(model1 <- fit_nonlinear(WSEw_obs[[k]][[m]], h=10))
      if (!is.null(model1)) # this does not take up much time, so no need to remove it
      {
        nl[[k]][[m]] <- model1
      } else
      {
        nl[[k]][[m]] <- NA
      }
    }
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
      if (!is.null(model1))
      {
        nlsb[[k]][[m]] <- model1
      } else
      {
        print("hi")
        nlsb[[k]][[m]] <- NA
      }
    }
  }
  nlsb_name <- paste0("nlsb/nlsb_", "r_", r, ".rds")
  saveRDS(nlsb, file = file.path(exp_dir, nlsb_name))
  return(0)
}