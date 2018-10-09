#' Linear Fit for WSE-w Relationship
#' 
#' Optionally uses Mersel's original method, which screens out "non-optimal" cross sections
#' @param WSEw WSEw data (at a given level of exposure)
#' @param mersel flag for Mersel method
#' @param thres threshold for Mersel method
#' @export

fit_linear <- function(WSEw,  mersel = FALSE, thres = NULL)
{

  if (!mersel)
  {
    lf <- lm(WSE~w, data = WSEw)
    return(lf)
    
  } else
  {
    tl <- test_linear(WSEw$WSE, WSEw$w, thres = thres)
    if (tl)
    {
      lf <- lm(WSE~w, data = WSEw)
      return(lf)
      
    } else
    {
      print("Not an optimal location")
      return(NULL)
    }
  }
}

# ------------------------------------------------------------------------------------------------------
#' Test Linear
#' 
#' Tests whether the linear method of Mersel et al. (2013) is appropriate for the WSE-w measurements
#' @export
test_linear <- function(WSE, w, thres = 0.015)
{
  # Tests if the linear method is appropriate for depth estimation
  maxdiff <- get_maxdiff(WSE,w)
  
  if (maxdiff < thres)
  {
    #print("Use the linear method")
    return(TRUE)
  } else
  {
    #print("Do not use the linear method")
    return(FALSE)
  }
  
}

