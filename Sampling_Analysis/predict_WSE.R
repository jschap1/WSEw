#' Predict WSE
#' 
#' Predicts WSE values that could be observed by SWOT based on a calibrated model
#' @details 
#' Logistic distribution as fit to Miss. Pool 21 Quincy gage, 1947-2018
#' @export
# THIS DID NOT WORK AT ALL. VERY BAD. 

predict_WSE <- function(N)
{
  
  # Logistic distribution parameters
  loc <- -3.597e-3
  scale <- 0.1199
  
  # Simulate data
  M <- 10000 # used for spin-up
  zp <- rlogis(N+M, loc, scale)
  
  # Un-difference
  z <- vector(length = N+M)
  z[1] <- 0
  for (t in 1:(N+M-1))
  {
    z[t+1] <- z[t] + zp[t+1]
  }
  # Remove the spin up
  z <- z[(M+1):(N+M)]
  
  # Re-season
  e <- vector(length = N) # residuals
  month.ind <- 1 # starts at January 1
  day.ind <- 1  
  for (i in 1:N)
  {
    e[i] <- z[i]*monthly.sd[month.ind]+monthly.mean[month.ind]
    
    day.ind <- day.ind + 1
    if (month.ind == 2 & day.ind > 28)
    {
      month.ind <- month.ind + 1 # ,6,9,11
    } else if ((month.ind == 4 & day.ind > 30) | 
               (month.ind == 6 & day.ind > 28) | 
               (month.ind == 9 & day.ind > 28)
               | (month.ind == 11 & day.ind > 28))
    {
      month.ind <- month.ind + 1
    } else if (day.ind > 31)
    {
      month.ind <- month.ind + 1
    }
    if (month.ind > 12) 
    {
      month.ind <- 1
    }
  }
  
  # Add the trend back in
  a <- coef(m1)[1]
  b <- coef(m1)[2]
  t2 <- (1:N) # not quite right, perhaps?
  hl <- e + a + b*t2
  H.pred <- as.numeric(exp(hl))
  
}