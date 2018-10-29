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