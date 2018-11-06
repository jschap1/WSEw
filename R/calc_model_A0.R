#' Calculate modeled flow area beneath the lowest observed WSE value
#' 
#' Assumes some model, such as linear, nonlinear, or slope break for the channel cross section
#' @export
#' @param model fitted model to the WSE-w data
#' @param type type of model. Default is linear. Options are linear, sb, sbm, nl, and nlsb
#' @param w1 minimum observed width value for the fit. Only required if fit is nonlinear
#' @param h1 minimum observed height value for the fit. Only required if fit is nonlinear
#' @details 
#' @return A0 cross-sectional area of flow beneath lowest observation
#' @example A0 <- calc_model_A0(lf1, type = "linear")

calc_model_A0 <- function(model, type, w1 = NULL, h1 = NULL, pos.only = FALSE)
{
  
  if (type == "linear")
  {
    
    w1 <- min(model$model$w)
    WSE0 <- min(model$model$WSE)
    z0 <- predict(model, newdata = data.frame(w=0))
    A <- 0.5*w1*(WSE0 - z0)
    
  }else if (type == "sb")
  {
    m1 <- model[[1]] # lower model
    w1 <- min(m1$model$w)
    WSE0 <- min(m1$model$WSE)
    z0 <- predict(m1, newdata = data.frame(w=0))
    A <- 0.5*w1*(WSE0 - z0)
    
  }else if (type == "sbm")
  {
    
    model <- model[[1]]
    w1 <- min(model$model$w)
    WSE0 <- min(model$model$WSE)
    z0 <- predict(model, newdata = data.frame(w=0))
    A <- 0.5*w1*(WSE0 - z0)
    
  }else if (type == "nl")
  {
    
    z0 <- model$z0
    a <- model$a
    s <- model$s
    A <- nl_area(w1, h1, z0, a, s)
    
  }else if (type == "nlsb")
  {
    
    z0 <- model$z0
    a1 <- model$a1
    s1 <- model$s1
    A <- nl_area(w1, h1, z0, a1, s1)

  }
  
  if (pos.only)
  {
    
    A[A<0] <- NA
    return(A)

  } else 
  {
    return(A)
  }

}

# --------------------------------------------------------------------------------------------------

#' Calculate area of a nonlinear cross section
#' @param s shape parameter for nonlinear fit (the exponent)
#' @param a multiplicative parameter for the nonlinear fit

nl_area <- function(wbf, WSEbf, z0, a, s)
{
  A <- (WSEbf - z0)*wbf - (a/(s+1))*wbf^(s+1)
  return(A)
}


