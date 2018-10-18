#' Predict parameters (parallel)
#'
#' Parallelized functions for main_d6 hydraulic parameter predictions in parallel with foreach
#' for WSE-w fits.
#' @export
#' @param r cross section number
#' @param exclude boolean that determines whether or not to exclude predictions from physically unrealistic model fits
#' @details Loads fitted model (linear, slope break, multiple slope break, nonlinear, or nonlinear slope break)
#' and uses it to predict hydraulic parameters
#' @return list of z0, A, W0, and A0 predictions
#' @example pred_vals <- pred_linear_par(r = 1)

# Linear
pred_linear_par <- function(r, exclude = FALSE)
{
  z0.l <- array(dim = c(n_exp_levels, M)) # initialize
  A.l <- array(dim = c(n_exp_levels, M))
  WP.l <- array(dim = c(n_exp_levels, M))
  A0.l <- array(dim = c(n_exp_levels, M))
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
        if (exclude == TRUE)
        {
          a <- coef(lf[[k]][[m]])[2] # slope coefficient
          if (a < 0)
          {
            next
          }
        }

        # z0.l[k,m] <- predict(lf[[k]][[m]], newdata = data.frame(w = 0))
        z0.l[k,m] <- as.numeric(coef(lf[[k]][[m]])[1])
        A.l[k,m] <- calc_model_A(lf[[k]][[m]], type = "linear")
        WP.l[k,m] <- calc_model_WP(lf[[k]][[m]], type = "linear")
        A0.l[k,m] <- calc_model_A0(lf[[k]][[m]], type = "linear", pos.only = FALSE)
      }
    }
  }
  return(list(z0 = z0.l, A = A.l, WP = WP.l, A0 = A0.l))
}


#' Slope break
#' 
#' Applies fit_slopebreak in parallel across multiple processors on one machine
#' @export
pred_sb_par <- function(r, exclude = FALSE)
{
  sb_name <- paste0("sb/sb_", "r_", r, ".rds")
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

        # Check that physical-realism constraints are satisfied
        # Must have non-negative slope
        if (exclude == TRUE)
        {
          a1 <- coef(sb[[k]][[m]][[1]])[2] # slope coefficient
          if (a1 < 0)
          {
            next
          }
        }

        # z0.sb[k,m] <- predict(sb[[k]][[m]][[1]], newdata = data.frame(w = 0))
        z0.sb[k,m] <- as.numeric(coef(sb[[k]][[m]][[1]])[1])
        A.sb[k,m] <- calc_model_A(sb[[k]][[m]], type = "sb")
        WP.sb[k,m] <- calc_model_WP(sb[[k]][[m]], type = "sb")
        A0.sb[k,m] <- calc_model_A0(sb[[k]][[m]], type = "sb", pos.only = FALSE)
      }
    }
  }
  return(list(z0 = z0.sb, A = A.sb, WP = WP.sb, A0 = A0.sb))
}

#' @export
#' Slope break (multiple)
#' 
#' Applies fit_slopebreak in parallel across multiple processors on one machine
pred_sbm_par <- function(r)
{
  sbm_name <- paste0("sbm/sbm_", "r_", r, ".rds")
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
        # try(z0.sbm[k,m] <- predict(sbm[[k]][[m]][[1]], newdata = data.frame(w = 0)))
        try(z0.sbm[k,m] <- as.numeric(coef(sbm[[k]][[m]][[1]])[1]))
        try(A.sbm[k,m] <- calc_model_A(sbm[[k]][[m]], type = "sbm"))
        try(WP.sbm[k,m] <- calc_model_WP(sbm[[k]][[m]], type = "sbm"))
        try(A0.sbm[k,m] <- calc_model_A0(sbm[[k]][[m]], type = "sbm", pos.only = TRUE))
      }
    }
  }
  return(list(z0 = z0.sbm, A = A.sbm, WP = WP.sbm, A0 = A0.sbm))
}


#' Nonlinear
#' 
#' Applies fit_nonlinear in parallel across multiple processors on one machine
#' @export
#' @param w1 minimum observed width, nr by n_exp_level by M array
#' @param h1 minimum observed height, nr by n_exp_level by M array

pred_nl_par <- function(r, WSEw, w1, h1, exclude = FALSE)
{
  nl_name <- paste0("nl/nl_", "r_", r, ".rds")
  if (file.exists(file.path(exp_dir, nl_name))) # error check in case no model was fit for this cross section
  {
    nl <- readRDS(file.path(exp_dir, nl_name))
  } else
  {
    return(list(z0 = NA, A = NA, WP = NA, A0 = NA))
  }
  z0.nl <- array(dim = c(n_exp_levels, M))
  A.nl <- array(dim = c(n_exp_levels, M))
  WP.nl <- array(dim = c(n_exp_levels, M))
  A0.nl <- array(dim = c(n_exp_levels, M))
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      if (class(nl[[k]][[m]]) == "nls")
      {

        # Check that physical-realism constraints are satisfied
        # Must have non-negative slope and cannot be concave down
        if (exclude == TRUE)
        {
          a <- as.numeric(coef(nl[[k]][[m]])[2]) # slope coefficient
          s <- as.numeric(coef(nl[[k]][[m]])[3]) # shape parameter
          cond1 <- (a>=0 & s>=1)
          cond2 <- (a<=0 & s<=1 & s>=0)
          if (!(cond1 | cond2))
          {
            next
          }
        }

        # z0.nl[k,m] <- predict(nl[[k]][[m]], newdata = data.frame(w = 0)) # this is problematic when s<0
        z0.nl[k,m] <- as.numeric(coef(nl[[k]][[m]])[1]) # this is better, will need to change it throughout
        A.nl[k,m] <- calc_model_A(nl[[k]][[m]], type = "nl", WSEw = WSEw[[r]])
        WP.nl[k,m] <- calc_model_WP(nl[[k]][[m]], type = "nl", w = WSEw[[r]]$w)
        A0.nl[k,m] <- calc_model_A0(nl[[k]][[m]], type = "nl", w1 = w1[r,k,m], h1 = h1[r,k,m], pos.only = FALSE)
      }
    }
  }
  return(list(z0 = z0.nl, A = A.nl, WP = WP.nl, A0 = A0.nl))
}

ww <- 0
6.664e+02 - 5.122e-04*ww^(-1.502e-01)


#' Nonlinear slope break
#' 
#' Applies fit_nlsb in parallel across multiple processors on one machine
#' @export
#' @param w1 minimum observed width, nr by n_exp_level by M array
#' @param h1 minimum observed height, nr by n_exp_level by M array
pred_nlsb_par <- function(r, WSEw, w1, h1, exclude = FALSE)
{
  nlsb_name <- paste0("nlsb/nlsb_", "r_", r, ".rds")
  if (file.exists(file.path(exp_dir, nlsb_name))) # error check in case no model was fit for this cross section
  {
    nlsb <- readRDS(file.path(exp_dir, nlsb_name))
  } else
  {
    return(list(z0 = NA, A = NA, WP = NA, A0 = NA))
  }
  z0.nlsb <- array(dim = c(n_exp_levels, M))
  A.nlsb <- array(dim = c(n_exp_levels, M))
  WP.nlsb <- array(dim = c(n_exp_levels, M))
  A0.nlsb <- array(dim = c(n_exp_levels, M))
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      if (class(nlsb[[k]][[m]][[1]]) == "nls")
      {

        # Check that physical-realism constraints are satisfied
        # Must have non-negative slope and cannot be concave down
        if (exclude == TRUE)
        {
          a1 <- as.numeric(coef(nlsb[[k]][[m]][[1]])[2]) # slope coefficient
          s1 <- as.numeric(coef(nlsb[[k]][[m]][[1]])[3]) # shape parameter
          cond1 <- (a1>=0 & s1>=1)
          cond2 <- (a1<=0 & s1<=1 & s1>=0)
          if (!(cond1 | cond2))
          {
            next
          }
        }

        # z0.nlsb[k,m] <- predict(nlsb[[k]][[m]][[1]], newdata = data.frame(w = 0))
        z0.nlsb[k,m] <- as.numeric(coef(nlsb[[k]][[m]][[1]])[1])
        A.nlsb[k,m] <- calc_model_A(nlsb[[k]][[m]], type = "nlsb", WSEw = WSEw[[r]]) # there may be a bug in the type = nlsb code here
        WP.nlsb[k,m] <- calc_model_WP(nlsb[[k]][[m]], type = "nlsb", w = WSEw[[r]]$w)
        A0.nlsb[k,m] <- calc_model_A0(nlsb[[k]][[m]], type = "nlsb", w1 = w1[r,k,m], h1 = h1[r,k,m], pos.only = FALSE)
      }
    }
  }
  return(list(z0 = z0.nlsb, A = A.nlsb, WP = WP.nlsb, A0 = A0.nlsb))
}


#' Clean by z0
#'
#' This removes the predictions with z0<h1
#' @export
#' @param pred list of predicted values
#' @param h1 array of minimum observed heights
clean_by_z0 <- function(pred, h1)
{
  for (r in 1:nr)
  {
    for (k in 1:n_exp_levels)
    {
      neglect_ind <- which(pred[[r]]$z0[k,] > h1[r,k,])
      pred[[r]]$z0[k,neglect_ind] <- NA
      pred[[r]]$A0[k,neglect_ind] <- NA
    }
  }
  return(pred)
}

#' Clean by A0
#' 
#' This removes the predictions with A0<0
#' @export
#' @details
clean_by_A0 <- function(pred, h1)
{
  for (r in 1:nr)
  {
    for (k in 1:n_exp_levels)
    {
      neglect_ind <- which(pred[[r]]$A0[k,] < 0)
      pred[[r]]$z0[k,neglect_ind] <- NA
      pred[[r]]$A0[k,neglect_ind] <- NA
    }
  }
  return(pred)
}
