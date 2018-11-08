#' Nonlinear Slope Break Fit
#' 
#' Fits a piecewise nonlinear curve with one breakpoint
#' There is one breakpoint, it is estimated using strucchange::breakpoints 
#' @param WSEw WSEw data (at a given level of exposure)
#' @param h minimum distance from either of the WSE-w curve where the breakpoint is constrained to be
#' @importFrom strucchange breakpoints
#' @importFrom minpack.lm nlsLM
#' @export
#' @examples
#' rWSEw <- readRDS("./Outputs/p4_newfits/rWSEw.rds")
#' r <- 5
#' WSEw <- observe(rWSEw[[r]], exposure = 0.4)
#' plot(WSE~w, WSEw)
#' nn <- length(WSEw$w)
#' fits <- fit_nlsb(WSEw)
#' sb.ind <- attributes(fits)$sb.ind
#' lines(WSEw$w[1:sb.ind], predict(fits[[1]]))
#' lines(WSEw$w[sb.ind:nn], predict(fits[[2]]))
#' Find an error case and see if it is handled better now
#' Did not manage to actually find such an error case
#' z0.error <- readRDS("./Outputs/p4_newfits/z0_error.rds")
#' z0.error.nlsb <- z0.error$nlsb
#' r <- 4
#' k <- 9
#' z0.error.nlsb[r,k,]
#' m <- 12
#' fit1 <- fit_nlsb(WSEw_obs[[k]][[m]])
#' @details See attributes(fit)$ef for error flags. 
#' Value 0 means no error, 
#' 1 means not enough data points, 
#' 2 means singular gradient at initial guess for nlsLM

fit_nlsb <- function(WSEw, h = 10, maxiter = 100)
{
  
  if (length(WSEw$WSE)<h) 
  {
    print("Not enough data points")
    fits <- NULL
    attributes(fits)$ef <- 1
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
      print("Using multiple initial guess method")
      
      ################################################
      
      init.guess <- make_init_guess(WSEw[1:i,])
      fits <- vector(length = dim(init.guess)[1], "list")
      SSE <- rep(NA, length = length(fits))
      
      # try performing each fit with the different sets of initial guesses
      for (p in 1:length(fits))
      {
        # if the algorithm fails to estimate parameters, return a warning message
        
        tryCatch({
          
          fits[[p]] <- nlsLM(WSE ~ a0 + a1*w^a2, 
                         start = c(a0 = init.guess$z0[p], a1 = init.guess$a[p], a2 = init.guess$s[p]), 
                         data = WSEw[1:i,], 
                         control = nls.lm.control(maxiter = maxiter)
          )
          
          
          SSE[p] <- sum(residuals(fits[[p]])^2)
          
        },
        
        error = function(e) 
        {
          print("error: nonlinear fit (multi) did not converge")
        }
        )
        
      }
      
      if (all(is.na(SSE)))
      {
        fits <- NULL
        attributes(fits)$ef <- 2
        return(fits)
      } else
      {
        fit.a <- fits[[which.min(SSE)]] # choose the best fit
      }
      
      ################################################
      
    } # end if !exists(fit.a) 
    
    if (!exists("fit.a"))
    {
      fits <- NULL
      attributes(fits)$ef <- 2
      return(fits)
    }
    
    if (exists("fit.a")) # if exists(fit.a)
    {
      # print("no problema with fit.a")
      
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
        print("Using multiple initial guess method")
        # use multiple initial guess method to try to get a result for fit.b
        
        ################################################
        
        init.guess <- make_init_guess(WSEw[1:i,], type = "nlsb", intercept = a0 + a1*xb^a2, xb = xb)
        fits <- vector(length = dim(init.guess)[1], "list")
        SSE <- rep(NA, length = length(fits))
        
        # try performing each fit with the different sets of initial guesses
        for (p in 1:length(fits))
        {
          # if the algorithm fails to estimate parameters, return a warning message
          
          tryCatch({
            
            fits[[p]] <- nlsLM(WSE ~ a0+a1*xb^a2 - b1*xb^b2 + b1*w^b2, 
                           start = c(b1 = init.guess$a[p], b2 = init.guess$s[p]), 
                           data = WSEw[(i):n,], 
                           control = nls.lm.control(maxiter = maxiter)
            ) 
            
            SSE[p] <- sum(residuals(fits[[p]])^2)
            
          },
          
          error = function(e) 
          {
            print("error: nonlinear fit (multi) did not converge")
          }
          )
          
        }
        
        if (all(is.na(SSE)))
        {
          fits <- NULL
          attributes(fits)$ef <- 2
          return(fits)
        } else
        {
          fit.b <- fits[[which.min(SSE)]] # choose the best fit
        }
        
        ################################################
        
      } # end if !exists(fit.b)
      
      if (!exists("fit.b"))
      {
        fits <- NULL
        attributes(fits)$ef <- 2
        return(fits)
      }
      
      if (exists("fit.b")) # if exists(fit.b)
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

  attributes(fits)$ef <- 0
  attributes(fits)$sb.ind <- sb.ind # adding the sb.ind as an output
  return(fits)
  
}
