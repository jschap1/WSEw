# --------------------------------------------------------------------------------------------------------------------
#' Calculate depth prediction error
#' 
#' Calculates depth prediction error, either absolute or as a fraction of bankfull depth
#' @param pred minimum bed elevation z0 predictions
#' @param truth true minimum bed elevation
#' @param h bankfull water surface elevation. For UMESC data, it is the reference WSE. Not needed for absolute error.
#' @param type can be "absolute" or "relative" or "mersel"
#' @export
#' @details Option "mersel" calculates standard error relative to bankfull depth, see Mersel et al. (2013)
#' @examples
#' refWSE <- 470*0.3048 # meters, for pool 21
#' z0.l.error.abs <- calc_depth_prediction_error(z0.l, z0.true.ra, type = "absolute")
#' z0.l.error.rel <- calc_depth_prediction_error(z0.l, z0.true.ra, h = refWSE, type = "relative")
#' load("./Outputs/Final/p21/z0_pred.rda")
#' z0.true.ra <- readRDS("./Outputs/Final/p21/z0_true_ra.rds")
#' truth <- z0.true.ra
#' pred <- z0.l
#' z0.L.SE <- calc_depth_prediction_error(z0.l, z0.true.ra, h = refWSE, type = "mersel")

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
    
  } else if (type == "mersel")
  {
    
    for (k in 1:n_exp_levels) 
    {
      SE <- array(dim = c(nr, n_exp_levels))
      for (r in 1:nr) # can vectorize for speed
      {
        for (k in 1:n_exp_levels)
        {
          SE[r,k] <- (truth[r] - mean(pred[r,k,], na.rm = TRUE))/(h - truth[r]) # SE of the mean prediction
          # SE[r,k] <- (truth[r] - median(pred[r,k,], na.rm = TRUE))/(h - truth[r]) # SE of the median prediction
          # e <- (truth[r] - pred[r,k,])/(h - truth[r])
          # e[is.na(e)] <- 0 # ignores NA values in the product
          # temp <- sqrt((1/M)*t(e)%*%e)
          # temp[temp==0] <- NA
          # SE[r,k] <- temp # rmse of the M replicates
        }
      }
      return(SE)
      
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
#' @param type can be "absolute" or "relative1" or "relative2"
#' @export
#' @examples A0.l.error.abs <- calc_A0_prediction_error(A0.l, A0.true.ra, type = "absolute")
#' A0.l.error.rel <- calc_A0_prediction_error(pred = A0.l, truth = A0.true.ra, A = A.true.ra, type = "relative1")
#' @details relative1 uses A for the denominator, but relative2 uses A0 for the denominator

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
    }
    
  } else if (type == "relative1")
  {
    for (k in 1:n_exp_levels)
    {
      pred.error[,k,] <- pred[,k,] - rep(truth[,k], M)
      for (r in 1:nr)
      {
        pred.error[r,k,] <- pred.error[r,k,]/(A[r])
      }
    }
  } else if (type == "relative2")
  {
    for (k in 1:n_exp_levels)
    {
      pred.error[,k,] <- pred[,k,] - rep(truth[,k], M)
      for (r in 1:nr)
      {
        pred.error[r,k,] <- pred.error[r,k,]/truth[r,k]
      }
    }
  }
  
  return(pred.error)
}

# --------------------------------------------------------------------------------------------------------------------
#' Calculate slope prediction error
#' 
#' Calculates slope prediction error, either absolute or _____
#' @param pred minimum bed elevation z0 predictions
#' @param truth true minimum bed elevation
#' @param type can be "absolute" or "relative"
#' @export
#' @examples s0.l.error.abs <- calc_A0_prediction_error(s0.l, s0.true.ra, type = "absolute")

calc_s0_prediction_error <- function(pred, truth, type)
{
  
  nr <- dim(pred)[1]
  n_exp_levels <- dim(pred)[2]
  M <- dim(pred)[3]
  
  pred.error <- array(dim = dim(pred))
  
  if (type == "absolute")
  {
    
    for (k in 1:n_exp_levels) 
    {
      pred.error[,k,] <- pred[,k,] - rep(truth, M)
    }
    
  } else if (type == "relative")
  {
    
    for (k in 1:n_exp_levels) 
    {
      pred.error[,k,] <- pred[,k,] - rep(truth, M)
      for (r in 1:nr)
      {
        pred.error[r,k,] <- pred.error[r,k,]/(truth[r])
      }
    }
    
  }
  return(pred.error)
}