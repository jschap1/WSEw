#' Nonlinear Fit for WSE-w Relationship
#' 
#' @param WSEw WSEw data (at a given level of exposure)
#' @param h minimum amount of data to perform fit
#' @importFrom minpack.lm nlsLM
#' @export

fit_nonlinear <- function(WSEw, h = 5)
{
  
  if (length(WSEw$WSE)<h) 
  {
    print("Not enough data points")
    return(NULL)
  }

  tryCatch(
    {
      fit <- nlsLM(WSE ~ b0 + b1*w^b2, start = list(b0 = min(WSEw$WSE), b1 = 1e-4, b2 = 2), data=WSEw)
      
      # fit <- nlsLM(WSE ~ b0 + b1*w^b2, start = list(b0 = min(WSEw_obs1$WSE), b1 = 1e-4, b2 = 2), data=WSEw_obs1)
      # 
      # init.guess <- make_init_guess(WSEw_obs1) # problem: as you s.init increases, so does the estimated parameter s
      # fit1 <- nlsLM(WSE ~ b0 + b1*w^b2, start = list(b0 = init.guess$z0[1], b1 = init.guess$a[1], b2 = init.guess$s[1]), data=WSEw_obs1)
      # fit2 <- nlsLM(WSE ~ b0 + b1*w^b2, start = list(b0 = init.guess$z0[2], b1 = init.guess$a[2], b2 = init.guess$s[2]), data=WSEw_obs1)
      # fit3 <- nlsLM(WSE ~ b0 + b1*w^b2, start = list(b0 = init.guess$z0[3], b1 = init.guess$a[3], b2 = init.guess$s[3]), data=WSEw_obs1)
      # fit4 <- nlsLM(WSE ~ b0 + b1*w^b2, start = list(b0 = init.guess$z0[4], b1 = init.guess$a[4], b2 = init.guess$s[4]), data=WSEw_obs1)
      # fit5 <- nlsLM(WSE ~ b0 + b1*w^b2, start = list(b0 = init.guess$z0[5], b1 = init.guess$a[5], b2 = init.guess$s[5]), data=WSEw_obs1)
      # fit6 <- nlsLM(WSE ~ b0 + b1*w^b2, start = list(b0 = init.guess$z0[6], b1 = init.guess$a[6], b2 = init.guess$s[6]), data=WSEw_obs1)
    
    }, 
    error = function(e) 
    {
      print("error: nonlinear fit did not converge")
      return(NULL)
    }
  )
  
  return(fit)

}
