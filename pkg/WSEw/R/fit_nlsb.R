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

fit_nlsb <- function(WSEw, h = 10)
{
  
  if (length(WSEw$WSE)<h) 
  {
    print("Not enough data points")
    return(NULL)
  }
  
  n <- dim(WSEw)[1] # number of data points

  LSE.best <- 1e6 # starting value
  for (i in (1+h):(n-h))
  {
    # for each possible breakpoint, find the LSE of the best fit
    
    # fit.a <- nlsLM(WSE ~ a0 + a1*w^a2, start = c(a0 = min(WSEw$WSE), a1 = 1e-4, a2 = 2), 
    #                data = WSEw[1:i,])
    
    tryCatch(
      {
        fit.a <- nlsLM(WSE ~ a0 + a1*w^a2, 
                       start = c(a0 = min(WSEw$WSE), a1 = 1e-4, a2 = 2), 
                       data = WSEw[1:i,])
      }, 
      error = function(e) 
      {
        print("error: nonlinear fit did not converge")
        fit.a <- NULL
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
        print("error: nonlinear fit did not converge")
        fit.b <- NULL
      }
    )
    
    # fit.b <- nlsLM(WSE ~ a0+a1*xb^a2 - b1*xb^b2 + b1*w^b2, 
    #                start = c(b1 = 1e-4, b2 = 2), 
    #                data = WSEw[(i):n,])
    
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
