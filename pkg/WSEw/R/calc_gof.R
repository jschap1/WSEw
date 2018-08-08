#' Calculate goodness of fit metrics
#' 
#' Calculates SSE, AIC, BIC, r2, mean absolute error, and z0 bias
#' for WSE-w fits.
#' @param model a model fit like the result of lm
#' @return SSE, AIC, BIC, r2, MAE (mean absolute error), mAE (maximum absolute error), z0.error

calc_gof <- function(model)
{
  
  SSE <- sum(model$residuals^2)
  AIC <- 2*2+21*log(SSE)+21*log(1/21) # check formula or just use built-in R formula AIC
  AIC <- AIC(model)
  BIC <- BIC(model)
  r2 <- summary(lf)$r.squared
  MAE <- mean(abs(lf$residuals))
  mAE <- max(abs(lf$residuals))
  
  z0.true <- lf$model$WSE[1]
  z0.pred <- predict(lf, newdata = data.frame(w = 0))
  z0.error <- z0.pred - z0.true
  
  gof <- data.frame(SSE = SSE, AIC = AIC, BIC = BIC, r2 = r2, 
                    MAE = MAE, mAE = mAE, z0.error = z0.error)
  return(gof)
}