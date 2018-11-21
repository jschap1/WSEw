#' Calculate modeled flow area beneath the lowest observed WSE value
#' 
#' Assumes some model, such as linear, nonlinear, or slope break for the channel cross section
#' @export
#' @param model fitted model to the WSE-w data
#' @param type type of model. Default is linear. Options are linear, sb, sbm, nl, and nlsb
#' @param WSEw_obs height and width observations data frame. Require for nl and nlsb models
#' @details Updated 11/21/2018 with a different method for choosing w1, h1 for A0 calculations. 
#' Chooses the lowest h1 > z0_pred and calculates w1 using the fitted model.
#' @return A0 cross-sectional area of flow beneath lowest observation
#' @example A0 <- calc_model_A0(lf1, type = "linear")

calc_model_A0 <- function(model, type, WSEw_obs = NULL, pos.only = FALSE)
{
  
  if (type == "linear")
  {
    
    a <- coef(model)[2]
    z0 <- coef(model)[1]
    h_obs <- model$model$WSE
    h1 <- min(h_obs[h_obs>z0])
    w1 <- (h1 - z0)/a
    A <- 0.5*w1*(h1 - z0)
    
  }else if (type == "sb")
  {
    
    m1 <- model[[1]] # lower model
    a <- coef(m1)[2]
    z0 <- coef(m1)[1]
    h_obs <- m1$model$WSE
    h1 <- min(h_obs[h_obs>z0])
    w1 <- (h1 - z0)/a
    A <- 0.5*w1*(h1 - z0)
  
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
    h_obs <- WSEw_obs$WSE
    h1 <- min(h_obs[h_obs>z0])
    w1 <- ((h1-z0)/a)^(1/s)
    A <- nl_area(w1, h1, z0, a, s)
    
  }else if (type == "nlsb")
  {
    
    z0 <- model$z0
    a1 <- model$a1
    s1 <- model$s1
    h_obs <- WSEw_obs$WSE
    h1 <- min(h_obs[h_obs>z0])
    w1 <- ((h1-z0)/a1)^(1/s1)
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


