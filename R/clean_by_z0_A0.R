#' Clean by z0
#'
#' This removes the predictions with z0<h1
#' @export
#' @param pred list of predicted values
#' @param h1 array of minimum observed heights
clean_by_z0 <- function(pred, h1)
{
  nr <- length(pred)
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
clean_by_A0 <- function(pred, h1)
{
  nr <- length(pred)
  for (r in 1:nr)
  {
    for (k in 1:n_exp_levels)
    {
      neglect_ind <- which(pred[[r]]$A0[k,] < 0)
      pred[[r]]$z0[k,neglect_ind] <- NA
      pred[[r]]$A0[k,neglect_ind] <- NA
      pred[[r]]$A[k,neglect_ind] <- NA
      pred[[r]]$WP[k,neglect_ind] <- NA
    }
  }
  return(pred)
}
