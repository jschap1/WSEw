# Backup

# Initialize
z0.l <- array(dim = c(nr, n_exp_levels, M))
z0.sb <- array(dim = c(nr, n_exp_levels, M))
z0.sbm <- array(dim = c(nr, n_exp_levels, M))
z0.nl <- array(dim = c(nr, n_exp_levels, M))
z0.nlsb <- array(dim = c(nr, n_exp_levels, M))

A.l <- array(dim = c(nr, n_exp_levels, M))
A.sb <- array(dim = c(nr, n_exp_levels, M))
A.sbm <- array(dim = c(nr, n_exp_levels, M))
A.nl <- array(dim = c(nr, n_exp_levels, M))
A.nlsb <- array(dim = c(nr, n_exp_levels, M))

WP.l <- array(dim = c(nr, n_exp_levels, M))
WP.sb <- array(dim = c(nr, n_exp_levels, M))
WP.sbm <- array(dim = c(nr, n_exp_levels, M))
WP.nl <- array(dim = c(nr, n_exp_levels, M))
WP.nlsb <- array(dim = c(nr, n_exp_levels, M))

A0.l <- array(dim = c(nr, n_exp_levels, M))
A0.sb <- array(dim = c(nr, n_exp_levels, M))
A0.sbm <- array(dim = c(nr, n_exp_levels, M))
A0.nl <- array(dim = c(nr, n_exp_levels, M))
A0.nlsb <- array(dim = c(nr, n_exp_levels, M))

# Linear
for (r in 1:nr) # takes 5 minutes for 50 cross sections at 19 exposure levels
{
  lf_name <- paste0("lf/lf_", "r_", r, ".rds")
  lf <- readRDS(file.path(exp_dir, lf_name))
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      # if statements handle cases where the model is NULL/no model was fit
      if (class(lf[[k]][[m]]) == "lm")
      {
        z0.l[r,k,m] <- predict(lf[[k]][[m]], newdata = data.frame(w = 0))
        A.l[r,k,m] <- calc_model_A(lf[[k]][[m]], type = "linear")
        WP.l[r,k,m] <- calc_model_WP(lf[[k]][[m]], type = "linear")
        A0.l[r,k,m] <- calc_model_A0(lf[[k]][[m]], type = "linear")
      }
    }
  }
  if (r %% 1 == 0)
  {
    print(paste("progress:", r, "of", nr))
  }
}

# Slope break
for (r in 1:nr)
{
  sb_name <- paste0("sb/sb_", "r_", r, ".rds")
  sb <- readRDS(file.path(exp_dir, sb_name))
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      if (class(sb[[k]][[m]][[1]]) == "lm")
      {
        z0.sb[r,k,m] <- predict(sb[[k]][[m]][[1]], newdata = data.frame(w = 0))
        A.sb[r,k,m] <- calc_model_A(sb[[k]][[m]], type = "sb")
        WP.sb[r,k,m] <- calc_model_WP(sb[[k]][[m]], type = "sb")
        A0.sb[r,k,m] <- calc_model_A0(sb[[k]][[m]], type = "sb")
      }
    }
  }
  if (r %% 1 == 0)
  {
    print(paste("progress:", r, "of", nr))
  }
}

# SBM
for (r in 1:nr)
{
  sbm_name <- paste0("sbm/sbm_", "r_", r, ".rds")
  sbm <- readRDS(file.path(exp_dir, sbm_name))
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      if (any(class(sbm[[k]][[m]][[1]])=="lm"))
      {
        z0.sbm[r,k,m] <- predict(sbm[[k]][[m]][[1]], newdata = data.frame(w = 0))
        A.sbm[r,k,m] <- calc_model_A(sbm[[k]][[m]], type = "sbm")
        WP.sbm[r,k,m] <- calc_model_WP(sbm[[k]][[m]], type = "sbm")
        A0.sbm[r,k,m] <- calc_model_A0(sbm[[k]][[m]], type = "sbm")
      }
    }
  }
  if (r %% 1 == 0)
  {
    print(paste("progress:", r, "of", nr))
  }
}

# Nonlinear
for (r in 1:nr)
{
  nl_name <- paste0("nl/nl_", "r_", r, ".rds")
  nl <- readRDS(file.path(exp_dir, nl_name))
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      if (class(nl[[k]][[m]]) == "nls")
      {
        z0.nl[r,k,m] <- predict(nl[[k]][[m]], newdata = data.frame(w = 0))
        A.nl[r,k,m] <- calc_model_A(nl[[k]][[m]], type = "nl", WSEw = rWSEw[[r]])
        WP.nl[r,k,m] <- calc_model_WP(nl[[k]][[m]], type = "nl", w = rWSEw[[r]]$w)
        A0.nl[r,k,m] <- calc_model_A0(nl[[k]][[m]], type = "nl", w0 = w0.ra[r,k])
      }
    }
  }
  if (r %% 1 == 0)
  {
    print(paste("progress:", r, "of", nr))
  }
}

# NLSB
for (r in 1:nr)
{
  nlsb_name <- paste0("nlsb/nlsb_", "r_", r, ".rds")
  nlsb <- readRDS(file.path(exp_dir, nlsb_name))
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      if (class(nlsb[[k]][[m]][[1]]) == "nls")
      {
        z0.nlsb[r,k,m] <- predict(nlsb[[k]][[m]][[1]], newdata = data.frame(w = 0))
        A.nlsb[r,k,m] <- calc_model_A(nlsb[[k]][[m]], type = "nlsb", WSEw = rWSEw[[r]]) # there may be a bug in the type = nlsb code here
        WP.nlsb[r,k,m] <- calc_model_WP(nlsb[[k]][[m]], type = "nlsb", w = rWSEw[[r]]$w)
        A0.nlsb[r,k,m] <- calc_model_A0(nlsb[[k]][[m]], type = "nlsb", w0 = w0.ra[r,k])
      }
    }
  }
  if (r %% 1 == 0)
  {
    print(paste("progress:", r, "of", nr))
  }
}

save(z0.l, z0.sb, z0.sbm, z0.nl, z0.nlsb, file = file.path(exp_dir, "z0_pred.rda"))
save(A.l, A.sb, A.sbm, A.nl, A.nlsb, file = file.path(exp_dir, "A_pred.rda"))
save(WP.l, WP.sb, WP.sbm, WP.nl, WP.nlsb, file = file.path(exp_dir, "WP_pred.rda"))
save(A0.l, A0.sb, A0.sbm, A0.nl, A0.nlsb, file = file.path(exp_dir, "A0_pred.rda"))
