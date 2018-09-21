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
    }, 
    error = function(e) 
    {
      print("error: nonlinear fit did not converge")
      return(NULL)
    }
  )
  
  return(fit)

}
