#' Nonlinear Slope Break Fit
#' 
#' Fits a piecewise nonlinear curve with one breakpoint
#' There is one breakpoint, it is estimated using strucchange::breakpoints 
#' @param WSEw WSEw data (at a given level of exposure)
#' @param h minimum distance from either of the WSE-w curve where the breakpoint is constrained to be
#' @importFrom strucchange breakpoints
#' @importFrom minpack.lm nlsLM
#' @export
#' @details Can be buggy. In particular, if the algorithm is unable to estimate model parameters, there can be an error where it says
#' that fit.b was not found. There should be a set of initial guesses that allow the solution to be found. 
#' The fix should just prevent fit.b from being called when it is not defined.

fit_nlsb <- function(WSEw, h = 10)
{
  
  if (length(WSEw$WSE)<h) 
  {
    print("Not enough data points")
    fits <- NULL
    return(fits)
  }
  
  n <- dim(WSEw)[1] # number of data points
  
  LSE.best <- 1e6 # starting value
  for (i in (1+h):(n-h))
  {
    # for each possible breakpoint, find the LSE of the best fit
    tryCatch(
      {
        fit.a <- nlsLM(WSE ~ a0 + a1*w^a2, 
                       start = c(a0 = min(WSEw$WSE), a1 = 1e-4, a2 = 2), 
                       data = WSEw[1:i,])
      }, 
      error = function(e) 
      {
        #print("error: nonlinear fit did not converge")
        fits <- NULL
        return(fits)
      }
    )
    
    a0 <- coef(fit.a)[1]
    a1 <- coef(fit.a)[2]
    a2 <- coef(fit.a)[3]
    xb <- WSEw$w[i]
    
    tryCatch(
      {
        fit.b <- nlsLM(WSE ~ a0+a1*xb^a2 - b1*xb^b2 + b1*w^b2, 
                       start = c(b1 = 1e-4, b2 = 2), 
                       data = WSEw[(i):n,])
      }, 
      error = function(e) 
      {
        #print("error: nonlinear fit did not converge")
        fits <- NULL
        return(fits)
      }
    )
    
    LSE <- sum(resid(fit.a)^2) + sum(resid(fit.b)^2)
    
    if (LSE<LSE.best)
    {
      # if LSE is low, record the estimate
      LSE.best <- LSE
      fits <- list(fit.a, fit.b)
      sb.ind <- i
    }
    
  }

  if (!exists("fits"))
  {
    print("fits does not exist. using fit_nlsb_try_multi")
    fits <- fit_nlsb_try_multi(WSEw, h = 10, maxiter = 100)
  } else if (is.null(fits))
  {
    print("fits is NULL. using fit_nlsb_try_multi")
    fits <- fit_nlsb_try_multi(WSEw, h = 10, maxiter = 100)
  }
  
  attributes(fits) <- list(sb.ind = sb.ind) # adding the sb.ind as an output
  return(fits)
  
}