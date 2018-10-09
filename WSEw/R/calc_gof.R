#' Calculate goodness of fit metrics
#'
#' Calculates SSE, AIC, BIC, r2, mean absolute error, and maximum absolute error
#' for WSE-w fits.
#' @export
#' @param model a model fit like the result of lm
#' @param type type of fit (if sb or nlsb, then the implementation has to change)
#' @details Calculates goodness-of-fit metrics for a model with respect to the fitted values.
#' R2 is not defined for nonlinear models, but it can be computed for a piecewise linear model
#' Piecewise AIC and BIC are computed by summing the AIC, BIC values of the original models. I am like 90% sure this is OK.
#' @return gof data frame containing SSE, AIC, BIC, r2, MAE (mean absolute error), and mAE (maximum absolute error)
#' @example gof <- calc_gof(lf1, type = "l)

calc_gof <- function(model, type)
{

  if (type == "l" | type == "sbm" | type == "nl")
  {

    if (type == "sbm")
    {
      model <- model[[1]]
    }

    SSE <- sum(residuals(model)^2)
    MAE <- mean(abs(residuals(model)))
    mAE <- max(abs(residuals(model)))

    AIC <- AIC(model)
    BIC <- BIC(model)

    if (type == "nl"){
      r2 <- NA # no r2 for a nonlinear model
    } else
    {
      r2 <- summary(model)$r.squared
    }

  } else if (type == "sb")
  {

    SSE <- sum(residuals(model[[1]])^2, residuals(model[[2]])^2)
    res1 <- abs(residuals(model[[1]]))
    res2 <- abs(residuals(model[[2]]))
    abs.residuals <- c(as.numeric(res1), as.numeric(res2))
    MAE <- mean(abs.residuals)
    mAE <- max(abs.residuals)

    # I am 90% sure AIC and BIC are linear operators
    AIC <- AIC(model[[1]]) + AIC(model[[2]])
    BIC <- BIC(model[[1]]) + BIC(model[[2]])

    # Calculate r2 for piecewise linear model
    SST1 <- SST(y = model[[1]]$model$WSE, y.bar = mean(model[[1]]$model$WSE))
    SST2 <- SST(y = model[[2]]$model$WSE, y.bar = mean(model[[2]]$model$WSE))
    r2_1 <- summary(model[[1]])$r.squared
    r2_2 <- summary(model[[2]])$r.squared
    r2 <- 1-(((1-r2_1)/(1+SST2/SST1))+((1-r2_2)/(1+SST1/SST2)))

  } else if (type == "nlsb")
  {

    SSE <- sum(residuals(model[[1]])^2, residuals(model[[2]])^2)
    res1 <- abs(residuals(model[[1]]))
    res2 <- abs(residuals(model[[2]]))
    abs.residuals <- c(as.numeric(res1), as.numeric(res2))
    MAE <- mean(abs.residuals)
    mAE <- max(abs.residuals)

    AIC <- AIC(model[[1]]) + AIC(model[[2]])
    BIC <- BIC(model[[1]]) + BIC(model[[2]])
    r2 <- NA # no r2 for a nonlinear model

  }

  gof <- data.frame(SSE = SSE, AIC = AIC, BIC = BIC, r2 = r2,
                    MAE = MAE, mAE = mAE)
  return(gof)
}

# ------------------------------------------------------------------------------------------

#' Batch goodness of fit metrics
#'
#' Use on a list of models
#' @param lf list of model objects
#' @return gof data frame containing goodness of fit measures for each model
#' @export

gof_batch <- function(model, type)
{
  nr <- length(model)
  for (r in 1:nr)
  {
    if (r==1)
    {
      gof <- calc_gof(model[[r]], type)
    } else
    {
      gof <- rbind(gof, calc_gof(model[[r]], type))
    }
    if (r%%5 == 0) # display progress
    {
      print(paste("Processed", r, "of", nr, "cross sections"))
    }
  }
  return(gof)
}

# ------------------------------------------------------------------------------------------

#' Calculate total sum of squares
#'
#' @param model model
#' @return SST
SST <- function(y, y.bar)
{
  SST <- 0
  n <- length(y)
  for (i in 1:n)
  {
    SST <- SST + (y[i]-y.bar)^2
  }
  return(SST)
}
