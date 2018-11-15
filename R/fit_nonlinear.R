#' Nonlinear Fit for WSE-w Relationship
#' 
#' @param WSEw WSEw data (at a given level of exposure)
#' @param h minimum amount of data to perform fit
#' @param maxiter maximum number of iterations for the parameter estimation algorithm
#' @importFrom minpack.lm nlsLM
#' @export
#' @details 
#' Uses a default initial guess, but if the fit is very bad, 
#' attempts several fits performed with different starting guesses 
#' and returns the best fit.
#' See attributes(fit)$ef for error flags. 
#' Value 0 means no error,
#' 1 means not enough data points, 
#' 2 means singular gradient at initial guess for nlsLM
#' @examples 
#' WSEw <- data.frame(WSE = NA, w = NA) # test case
#' model1 <- fit_nonlinear(WSEw, h = 5, maxiter = 100)

fit_nonlinear <- function(WSEw, h = 5, maxiter = 100)
{
  
  try_multi <- FALSE # flag for using the multiple initial guess method
  
  if (length(WSEw$WSE)<h) 
  {
    print("Not enough data points")
    fit <- NULL
    attributes(fit)$ef <- 1
    return(fit)
  }

  # General initial guess assuming s = 2
  # Try this method. If it fails, switch to the multiple initial guess method.
  lf <- lm(WSE ~ I(w^2), data = WSEw)
  a.init <- as.numeric(coef(lf)[2])
  init.guess <- list(z0 = min(WSEw$WSE), 
                     a = a.init, 
                     s = 2)
  tryCatch(
    {
      fit <- nlsLM(WSE ~ b0 + b1*w^b2, 
                   start = list(b0 = init.guess$z0, b1 = init.guess$a, b2 = init.guess$s), 
                   data=WSEw, 
                   control = nls.lm.control(maxiter = maxiter))
    }, 
    error = function(e) 
    {
      print("error: nonlinear fit did not converge")
    }
  )
  
  if (!exists("fit"))
  {
    try_multi <- TRUE
  } else
  {
    attributes(fit)$ef <- 0
  }
  
  if (try_multi)
  {
    
    init.guess <- make_init_guess(WSEw)
    fits <- vector(length = dim(init.guess)[1], "list")
    SSE <- rep(NA, length = length(fits))
    
    # try performing each fit with the different sets of initial guesses
    for (i in 1:length(fits))
    {
      # if the algorithm fails to estimate parameters, return a warning message
      
      tryCatch({
        
        fits[[i]] <- nlsLM(WSE ~ b0 + b1*w^b2, 
                           start = list(b0 = init.guess$z0[i], b1 = init.guess$a[i], b2 = init.guess$s[i]), 
                           data=WSEw, 
                           control = nls.lm.control(maxiter = maxiter)) 
        SSE[i] <- sum(residuals(fits[[i]])^2)
        
      },
      
      error = function(e) 
      {
        print("error: nonlinear fit (multi) did not converge")
      }
      )
      
    }
    
    if (all(is.na(SSE)))
    {
      fit <- NULL
      attributes(fit)$ef <- 2
    } else
    {
      fit <- fits[[which.min(SSE)]] # choose the best fit
      attributes(fit)$ef <- 0
    }

  }
  
  return(fit)

}
