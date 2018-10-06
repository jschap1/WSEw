# --------------------------------------------------------------------------------------------------------------------
#' Calculate prediction error
#' 
#' Calculates hydraulic parameter prediction error
#' @param truth
#' @param predicted
#' @param A0 flag if using A0
#' @export

calc_prediction_error <- function(pred, truth, A0 = FALSE)
{
  
  # nr <- dim(pred)[1]
  n_exp_levels <- dim(pred)[2]
  M <- dim(pred)[3]
  
  pred.error <- array(dim = dim(pred))
  
  if (!A0)
  {
    
    for (k in 1:n_exp_levels)
    {
      pred.error[,k,] <- pred[,k,] - rep(truth, M)
    }

  } else
  {
    for (k in 1:n_exp_levels)
    {
      pred.error[,k,] <- pred[,k,] - rep(truth[,k], M)
    }
  }
  
  return(pred.error)
  
}

