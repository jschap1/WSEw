# --------------------------------------------------------------------------------------------------------------------
#' Calculate depth prediction error
#' 
#' Calculates depth prediction error, either absolute or as a fraction of bankfull depth
#' @param pred minimum bed elevation z0 predictions
#' @param truth true minimum bed elevation
#' @param h bankfull water surface elevation. For UMESC data, it is the reference WSE. Needed for relative error only.
#' @param type can be "absolute" or "relative"
#' @export
#' @examples
#' refWSE <- 470*0.3048 # meters, for pool 21
#' z0.l.error.abs <- calc_depth_prediction_error(z0.l, z0.true.ra, type = "absolute")
#' z0.l.error.rel <- calc_depth_prediction_error(z0.l, z0.true.ra, h = refWSE, type = "relative")

calc_depth_prediction_error <- function(pred, truth, h = NULL, type)
{
  
  nr <- dim(pred)[1]
  n_exp_levels <- dim(pred)[2]
  M <- dim(pred)[3]
  
  pred.error <- array(dim = dim(pred))
  
  if (type == "absolute")
  {
    
    for (k in 1:n_exp_levels) 
    {
      pred.error[,k,] <- rep(truth, M) - pred[,k,]
    }
    
  } else if (type == "relative")
  {
    
    for (k in 1:n_exp_levels) 
    {
      pred.error[,k,] <- rep(truth, M) - pred[,k,]
      for (r in 1:nr)
      {
        pred.error[r,k,] <- pred.error[r,k,]/(h - truth[r])
      }
    }
    
  }
  return(pred.error)
}
  
# --------------------------------------------------------------------------------------------------------------------
#' Calculate submerged area prediction error
#' 
#' Calculates submerged area prediction error, either absolute or as a fraction of bankfull flow area
#' @param pred minimum bed elevation z0 predictions
#' @param truth true minimum bed elevation
#' @param A bankfull flow area
#' @param type can be "absolute" or "relative"
#' @export
#' @examples A0.l.error.abs <- calc_A0_prediction_error(A0.l, A0.true.ra, type = "absolute")
#' A0.l.error.rel <- calc_A0_prediction_error(pred = A0.l, truth = A0.true.ra, A = A.true.ra, type = "relative")

calc_A0_prediction_error <- function(pred, truth, A = NULL, type)
{
  
  nr <- dim(pred)[1]
  n_exp_levels <- dim(pred)[2]
  M <- dim(pred)[3]
  
  pred.error <- array(dim = dim(pred))
  
  if (type == "absolute")
  {
    for (k in 1:n_exp_levels)
    {
      pred.error[,k,] <- pred[,k,] - rep(truth[,k], M)
      # for (r in 1:nr)
      # {
      #   pred.error[r,k,] <- pred.error[r,k,]/truth[r,k]
      # }
    }
    
  } else if (type == "relative")
  {
    for (k in 1:n_exp_levels)
    {
      pred.error[,k,] <- pred[,k,] - rep(truth[,k], M)
      for (r in 1:nr)
      {
        pred.error[r,k,] <- pred.error[r,k,]/(A[r])
      }
    }
  }
  
  return(pred.error)
}