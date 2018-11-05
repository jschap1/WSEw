#' Nonlinear Fit for WSE-w Relationship
#' 
#' @param WSEw WSEw data (at a given level of exposure)
#' @param h minimum amount of data to perform fit
#' @importFrom minpack.lm nlsLM
#' @export
#' @details Returns the best fit, out of several fits performed with different starting guesses

fit_nonlinear <- function(WSEw, h = 5, maxiter = 100)
{
  
  if (length(WSEw$WSE)<h) 
  {
    print("Not enough data points")
    return(NULL)
  }

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
      print("error: nonlinear fit did not converge")
    }
    )
    
  }
  
  if (all(is.na(SSE)))
  {
    fit <- NULL
  } else
  {
    fit <- fits[[which.min(SSE)]] # choose the best fit
  }
  
  return(fit)

}
