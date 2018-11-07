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
#' @examples
#' WSEw <- rWSEw[[3]]
#' fits <- fit_nlsb(WSEw)

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
  hmin <- h/2 # minimum distance the breakpoint is assumed to be from either side
  for (i in (1+hmin):(n-hmin))
  {
    # for each possible breakpoint, find the LSE of the best fit
    try(
      {
        # General initial guess assuming s = 2
        # Try this method. If it fails, switch to the multiple initial guess method.
        lf <- lm(WSE ~ I(w^2), data = WSEw[1:i,])
        a.init <- as.numeric(coef(lf)[2])
        init.guess <- list(z0 = min(WSEw[1:i,]$WSE), 
                           a = a.init, 
                           s = 2)
        
        fit.a <- nlsLM(WSE ~ a0 + a1*w^a2, 
                       start = c(a0 = init.guess$z0, a1 = init.guess$a, a2 = init.guess$s), 
                       data = WSEw[1:i,], 
                       control = nls.lm.control(maxiter = maxiter)
                       )
      }, 
    ) # end try
    
    if (!exists("fit.a"))
    {
      
      # if the first fit failed, try the multiple initial guess method
      print("fit.a does not exist")
      
    } # end if !exists(fit.a) 
    else # if exists(fit.a)
    {
      a0 <- as.numeric(coef(fit.a)[1])
      a1 <- as.numeric(coef(fit.a)[2])
      a2 <- as.numeric(coef(fit.a)[3])
      xb <- WSEw$w[i]
      
      try(
        {
          intercept <- a0 + a1*xb^a2
          lf <- lm(I(WSE - intercept) ~ 0 + I(w^2 - xb^2), data = WSEw[(i):n,])
          a.init <- as.numeric(coef(lf)[1])
          init.guess <- list(a = a.init, s = 2)
          fit.b <- nlsLM(WSE ~ a0+a1*xb^a2 - b1*xb^b2 + b1*w^b2, 
                         start = c(b1 = init.guess$a, b2 = init.guess$s), 
                         data = WSEw[(i):n,], 
                         control = nls.lm.control(maxiter = maxiter)
                         )
        }, 
      ) # end try
      
      if(!exists("fit.b"))
      {
        print("fit.b does not exist")
        
      } # end if !exists(fit.b)
      else # if exists(fit.b)
      {
        
        LSE <- sum(resid(fit.a)^2) + sum(resid(fit.b)^2)
        
        if (LSE<LSE.best)
        {
          # if LSE is low, record the estimate
          LSE.best <- LSE
          fits <- list(fit.a, fit.b)
          sb.ind <- i
        }
        
      } # end else if exists(fit.b)

    } # end else if exists(fit.a)
    
  } # end for loop

  attributes(fits) <- list(sb.ind = sb.ind) # adding the sb.ind as an output
  return(fits)
  
}
