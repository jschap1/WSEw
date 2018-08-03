#' Nonlinear Fit for WSE-w Relationship
#' 
#' @param WSEw WSEw data (at a given level of exposure)
#' @export

fit_nonlinear <- function(WSEw)
{
  
  if (length(WSEw$WSE)<2*5) 
  {
    print("Not enough data points")
    return(NULL)
  }

  tryCatch(
    {

      fit <- nlsLM(WSE ~ b0 + b1*w^b2, start = list(b0 = min(WSEw$WSE), b1 = 0.01, b2 = 1), data=WSEw)
      
    }, 
    error = function(e) 
    {
      print("error: nonlinear fit did not converge")
      return(NULL)
    }
  )
  
  return(fit)

}
