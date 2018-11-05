#' Nonlinear Slope Break Fit
#' 
#' Fits a piecewise nonlinear curve with one breakpoint
#' There is one breakpoint, it is estimated using strucchange::breakpoints 
#' @param WSEw WSEw data (at a given level of exposure)
#' @param h minimum distance from either of the WSE-w curve where the breakpoint is constrained to be
#' @importFrom strucchange breakpoints
#' @importFrom minpack.lm nlsLM
#' @export
#' @details Eventually, it would be good to improve the initial guesses for the parameters, say, using linearization as described here:
#' https://stats.stackexchange.com/questions/160552/why-is-nls-giving-me-singular-gradient-matrix-at-initial-parameter-estimates
#' Can be buggy. In particular, if the algorithm is unable to estimate model parameters, there can be an error where it says
#' that fit.b was not found. There should be a set of initial guesses that allow the solution to be found. 
#' The fix should just prevent fit.b from being called when it is not defined.

fit_nlsb <- function(WSEw, h = 10, maxiter = 100)
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
    
    init.guess <- make_init_guess(WSEw[1:i,])
    fits <- vector(length = dim(init.guess)[1], "list")
    SSE <- rep(NA, length = length(fits))
    
    # try performing each fit with the different sets of initial guesses
    for (i in 1:length(fits))
    {
      tryCatch({
        fits[[i]] <- nlsLM(WSE ~ a0 + a1*w^a2, 
                           start = list(a0 = init.guess$z0[i], a1 = init.guess$a[i], a2 = init.guess$s[i]), 
                           data=WSEw[1:i,],
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
      fit.a <- NULL
    } else
    {
      fit.a <- fits[[which.min(SSE)]] # choose the best fit
    }
    
    a0 <- coef(fit.a)[1]
    a1 <- coef(fit.a)[2]
    a2 <- coef(fit.a)[3]
    xb <- WSEw$w[i]
    
    # ------------------------------------------------------------------------------------------------
    
    init.guess <- make_init_guess(WSEw[(i):n,])
    fits <- vector(length = dim(init.guess)[1], "list")
    SSE <- rep(NA, length = length(fits))
    
    # try performing each fit with the different sets of initial guesses
    for (i in 1:length(fits))
    {
      tryCatch({
        fits[[i]] <- nlsLM(WSE ~ a0+a1*xb^a2 - b1*xb^b2 + b1*w^b2, 
                           start = list(b1 = init.guess$a[i], b2 = init.guess$s[i]), 
                           data=WSEw[(i):n,],
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
      fit.b <- NULL
    } else
    {
      fit.b <- fits[[which.min(SSE)]] # choose the best fit
    }
    
    LSE <- sum(resid(fit.a)^2) + sum(resid(fit.b)^2)
    
    if (LSE<LSE.best)
    {
      # if LSE is low, record the estimate
      LSE.best <- LSE
      fits <- list(fit.a, fit.b)
      sb.ind <- i
    }
    
  }
  
  attributes(fits) <- list(sb.ind = sb.ind) # adding the sb.ind as an output
  return(fits)
  
}
