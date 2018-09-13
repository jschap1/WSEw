# Parallelized fitting functions for main_d6 model fitting in parallel with foreach

# Fit linear model
fit_linear_par <- function(r)
{
  lf <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels)
  {
    lf[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
      model1 <- fit_linear(WSEw_obs[[r]][[k]][[m]])
      if (!is.null(model1))
      {
        lf[[k]][[m]] <- model1
      } else
      {
        lf[[k]][[m]] <- model1 <- NA
      }
    }
  }
  lf_name <- paste0("lf/lf_", "r_", r, ".rds")
  saveRDS(lf, file = file.path(exp_dir, lf_name))
  return(0)
}

# Fit SB model
fit_sb_par <- function(r)
{
  sb <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels)
  {
    sb[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
      try(model1 <- fit_slopebreak(WSEw_obs[[r]][[k]][[m]], multiple_breaks = FALSE, continuity = TRUE)) # just in case it decides to throw an error
      if (!is.null(model1))
      {
        sb[[k]][[m]] <- model1
      } else
      {
        sb[[k]][[m]] <- NA
      }
    }
  }
  sb_name <- paste0("sb/sb_", "r_", r, ".rds")
  saveRDS(sb, file = file.path(exp_dir, sb_name))
  return(0)
}

# To do: use the advice here: https://stackoverflow.com/questions/12135400/errors-in-segmented-package-breakpoints-confusion
# This will likely allow sbm fits to work more often, by restarting multiple times

# Fit SBM model
fit_sbm_par <- function(r)
{
  sbm <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels)
  {
    sbm[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
      try(model1 <- fit_slopebreak(WSEw_obs[[r]][[k]][[m]], multiple_breaks = TRUE, continuity = TRUE)) # sometimes this throws errors
      if (!is.null(model1))
      {
        sbm[[k]][[m]] <- model1
      } else
      {
        sbm[[k]][[m]] <- NA
      }
    }
  }
  sbm_name <- paste0("sbm/sbm_", "r_", r, ".rds")
  saveRDS(sbm, file = file.path(exp_dir, sbm_name))
  return(0)
}

# Fit nonlinear model
fit_nl_par <- function(r)
{
  nl <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels)
  {
    nl[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
      try(model1 <- fit_nonlinear(WSEw_obs[[r]][[k]][[m]]))
      if (!is.null(model1))
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

# Fit NLSB model
fit_nlsb_par <- function(r)
{
  nlsb <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels)
  {
    nlsb[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
      try(model1 <- fit_nlsb(WSEw_obs[[r]][[k]][[m]]))
      if (!is.null(model1))
      {
        nlsb[[k]][[m]] <- model1
      } else
      {
        nlsb[[k]][[m]] <- NA
      }
    }
  }
  nlsb_name <- paste0("nlsb/nlsb_", "r_", r, ".rds")
  saveRDS(nlsb, file = file.path(exp_dir, nlsb_name))
  return(0)
}