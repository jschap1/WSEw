# Parallelized functions for main_d6 hydraulic parameter predictions in parallel with foreach

# Linear
pred_linear_par <- function(r)
{
  z0.l <- array(dim = c(n_exp_levels, M)) # initialize
  A.l <- array(dim = c(n_exp_levels, M))
  WP.l <- array(dim = c(n_exp_levels, M))
  A0.l <- array(dim = c(n_exp_levels, M))
  lf_name <- paste0("lf/lf_", "r_", r, "_test.rds") # load fitted model
  lf <- readRDS(file.path(exp_dir, lf_name))
  for (k in 1:n_exp_levels)
  { 
    for (m in 1:M)
    {
      # if statements handle cases where the model is NULL/no model was fit
      if (class(lf[[k]][[m]]) == "lm")
      {
        z0.l[k,m] <- predict(lf[[k]][[m]], newdata = data.frame(w = 0))
        A.l[k,m] <- calc_model_A(lf[[k]][[m]], type = "linear")
        WP.l[k,m] <- calc_model_WP(lf[[k]][[m]], type = "linear")
        A0.l[k,m] <- calc_model_A0(lf[[k]][[m]], type = "linear")
      }
    }
  }
  return(list(z0 = z0.l, A = A.l, WP = WP.l, A0 = A0.l))
}

# Slope break
pred_sb_par <- function(r)
{
  sb_name <- paste0("sb/sb_", "r_", r, "_test.rds")
  sb <- readRDS(file.path(exp_dir, sb_name))
  z0.sb <- array(dim = c(n_exp_levels, M))
  A.sb <- array(dim = c(n_exp_levels, M))
  WP.sb <- array(dim = c(n_exp_levels, M))
  A0.sb <- array(dim = c(n_exp_levels, M))
  for (k in 1:n_exp_levels)
  { 
    for (m in 1:M)
    {
      if (class(sb[[k]][[m]][[1]]) == "lm")
      {
        z0.sb[k,m] <- predict(sb[[k]][[m]][[1]], newdata = data.frame(w = 0))
        A.sb[k,m] <- calc_model_A(sb[[k]][[m]], type = "sb")
        WP.sb[k,m] <- calc_model_WP(sb[[k]][[m]], type = "sb")
        A0.sb[k,m] <- calc_model_A0(sb[[k]][[m]], type = "sb")
      }
    }
  }
  return(list(z0 = z0.sb, A = A.sb, WP = WP.sb, A0 = A0.sb))
}

# SBM
pred_sbm_par <- function(r)
{
  sbm_name <- paste0("sbm/sbm_", "r_", r, "_test.rds")
  print(sbm_name)
  if (file.exists(file.path(exp_dir, sbm_name))) # error check in case no model was fit for this cross section
  {
    sbm <- readRDS(file.path(exp_dir, sbm_name))
  } else
  {
    return(list(z0 = NA, A = NA, WP = NA, A0 = NA))
  }
  z0.sbm <- array(dim = c(n_exp_levels, M))
  A.sbm <- array(dim = c(n_exp_levels, M))
  WP.sbm <- array(dim = c(n_exp_levels, M))
  A0.sbm <- array(dim = c(n_exp_levels, M))  
  for (k in 1:n_exp_levels)
  { 
    for (m in 1:M)
    {
      if (any(class(sbm[[k]][[m]][[1]])=="lm"))
      {
        try(z0.sbm[k,m] <- predict(sbm[[k]][[m]][[1]], newdata = data.frame(w = 0)))
        try(A.sbm[k,m] <- calc_model_A(sbm[[k]][[m]], type = "sbm"))
        try(WP.sbm[k,m] <- calc_model_WP(sbm[[k]][[m]], type = "sbm"))
        try(A0.sbm[k,m] <- calc_model_A0(sbm[[k]][[m]], type = "sbm"))
      }
    }
  }
  return(list(z0 = z0.sbm, A = A.sbm, WP = WP.sbm, A0 = A0.sbm))
}

# Nonlinear
pred_nl_par <- function(r)
{
  z0.nl <- array(dim = c(n_exp_levels, M))
  A.nl <- array(dim = c(n_exp_levels, M))
  WP.nl <- array(dim = c(n_exp_levels, M))
  A0.nl <- array(dim = c(n_exp_levels, M))  
  nl_name <- paste0("nl/nl_", "r_", r, "_test.rds")
  nl <- readRDS(file.path(exp_dir, nl_name))
  for (k in 1:n_exp_levels)
  { 
    for (m in 1:M)
    {
      if (class(nl[[k]][[m]]) == "nls")
      {
        z0.nl[k,m] <- predict(nl[[k]][[m]], newdata = data.frame(w = 0))
        A.nl[k,m] <- calc_model_A(nl[[k]][[m]], type = "nl", WSEw = xWSEw[[r]])
        WP.nl[k,m] <- calc_model_WP(nl[[k]][[m]], type = "nl", w = xWSEw[[r]]$w)
        A0.nl[k,m] <- calc_model_A0(nl[[k]][[m]], type = "nl", w0 = w0.ra[r,k])
      }
    }
  }
  return(list(z0 = z0.nl, A = A.nl, WP = WP.nl, A0 = A0.nl))
}

# NLSB
pred_nlsb_par <- function(r)
{
  z0.nlsb <- array(dim = c(n_exp_levels, M))
  A.nlsb <- array(dim = c(n_exp_levels, M))
  WP.nlsb <- array(dim = c(n_exp_levels, M))
  A0.nlsb <- array(dim = c(n_exp_levels, M))  
  nlsb_name <- paste0("nlsb/nlsb_", "r_", r, "_test.rds")
  nlsb <- readRDS(file.path(exp_dir, nlsb_name))
  for (k in 1:n_exp_levels)
  { 
    for (m in 1:M)
    {
      if (class(nlsb[[k]][[m]][[1]]) == "nls")
      {
        z0.nlsb[k,m] <- predict(nlsb[[k]][[m]][[1]], newdata = data.frame(w = 0))
        A.nlsb[k,m] <- calc_model_A(nlsb[[k]][[m]], type = "nlsb", WSEw = xWSEw[[r]]) # there may be a bug in the type = nlsb code here
        WP.nlsb[k,m] <- calc_model_WP(nlsb[[k]][[m]], type = "nlsb", w = xWSEw[[r]]$w)
        A0.nlsb[k,m] <- calc_model_A0(nlsb[[k]][[m]], type = "nlsb", w0 = w0.ra[r,k])
      }
    }
  }
  return(list(z0 = z0.nlsb, A = A.nlsb, WP = WP.nlsb, A0 = A0.nlsb))
}



