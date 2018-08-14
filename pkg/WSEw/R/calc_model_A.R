#' Calculate modeled bankfull flow area
#' 
#' Assumes some model, such as linear, nonlinear, or slope break for the channel cross section
#' @param model fitted model to the WSE-w data
#' @param type type of model. Default is linear. Options are linear, sb, sbm, nl, and nlsb
#' @param WSEw WSE-w data for the fit. Only required if fit is nonlinear
#' @export
#' @return A cross-sectional area of flow
#' @examples 
#' A.l <- calc_modelA(lf[[1]], type = "linear")
#' A.sb <- calc_modelA(sbm[[1]], type = "sbm")
#' @details 

calc_model_A <- function(model, type, WSEw = NULL)
{
  
  if (type == "linear")
  {
    
    wbf <- max(model$model$w)
    WSEbf <- predict(model, newdata = data.frame(w=wbf))
    z0 <- predict(model, newdata = data.frame(w=0))
    A <- 0.5*wbf*(WSEbf - z0)
    
  }else if (type == "sb") # it would be more straightforward to use segmented for everything, including the single SB model
  {
    m1 <- model[[1]] # lower model
    m2 <- model[[2]] # upper model 
    wsb <- max(m1$model$w)
    wbf <- max(m2$model$`I(w - wb)`) + wsb # accounting for the offset
    WSEsb <- max(fitted(m1))
    WSEbf <- max(fitted(m2))
    z0 <- predict(m1, newdata = data.frame(w=0))
    A <- sb_area(wsb, WSEsb, z0, WSEbf, wbf)   
    
  }else if (type == "sbm")
  {
    
    model <- model[[1]]
    
    wbf <- max(model$model$w)
    # WSEbf <- max(model$model$WSE)
    WSEbf <- max(fitted(model))
    
    w.brk <- as.numeric(summary(model)$psi[,2]) # width breakpoints, in increasing order
    WSE.brk <- predict(model, newdata = data.frame(w = w.brk))
    
    w <- c(w.brk, wbf)
    WSE <- c(WSE.brk, WSEbf)
    # z0 <- min(model$model$WSE)
    z0 <- predict(model, newdata = data.frame(w = 0))
    
    A <- sbm_area(w, WSE, z0) 
    
  }else if (type == "nl")
  {
    
    wbf <- max(WSEw$w)
    WSEbf <- max(fitted(model))
    z0 <- predict(model, newdata = data.frame(w=0))
    a <- as.numeric(coef(model)[2])
    s <- as.numeric(coef(model)[3])
    A <- nl_area(wbf, WSEbf, z0, a, s)
    
  }else if (type == "nlsb")
  {
    
    m1 <- model[[1]] # lower model
    m2 <- model[[2]] # upper model 
    
    wbf <- max(WSEw$w)
    WSEbf <- max(fitted(m2))
    z0 <- predict(m1, newdata = data.frame(w=0))
    
    sb.ind <- as.numeric(attributes(model)) # index of slope break
    wsb <- WSEw$w[sb.ind] # width at slope break

    a1 <- as.numeric(coef(m1)[2])
    a2 <- as.numeric(coef(m2)[1])
    s1 <- as.numeric(coef(m1)[3])
    s2 <- as.numeric(coef(m2)[2])
    a <- c(a1,a2)
    s <- c(s1,s2)

    A <- nlsb_area3(wbf, WSEbf, z0, a, s, wsb)
  }
  
  return(A)
}


# --------------------------------------------------------------------------------------------------

#' Calculates area of a slope break cross section
 
sb_area <- function(wsb, WSEsb, z0, WSEbf, wbf)   
{
  A <- 0.5*(WSEsb-z0)*wsb + 0.5*(wsb+wbf)*(WSEbf - WSEsb)
  return(A)
}

# --------------------------------------------------------------------------------------------------

#' Calculates area of a multi-slope break cross section
#' @param w width
#' @param WSE water surface elevation. The lowest value is z0 and the highest is bankfull WSE.

sbm_area <- function(w, WSE, z0)
{
  n <- length(w) # number of segments in the model
  sum1 <- 0
  for (i in 1:(n-1))
  {
    sum1 <- sum1 + (w[i]+w[i+1])*(WSE[i+1]-WSE[i])
  }
  A <- 0.5*(WSE[1]-z0)*w[1] + 0.5*sum1
  return(A)
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

# --------------------------------------------------------------------------------------------------

#' Calculate area of a nonlinear slope break cross section
#' @param s shape parameters for nonlinear fits (the exponents)
#' @param a multiplicative parameters for the nonlinear fits

# nlsb_area <- function(wbf, WSEbf, z0, a, s, wsb)
# {
#   # Separating by term to keep the code looking clean
#   t1 <- (WSEbf - z0)*wbf # term 1
#   t2 <- a[1]*s[1]*wsb^(s[1]+1)/(s[1]+1) # term 2
#   t3 <- -a[1]*wbf*wsb^(s[1])
#   t4 <- (-a[2]/(s[2]+1))*(wbf^(s[2]+1) - wsb^(s[2]+1))
#   A <- sum(t1, t2, t3, t4)
#   return(A)
# }
# 
# nlsb_area2 <- function(wbf, WSEbf, z0, a, s, wsb)
# {
#   t1 <- WSEbf*wbf
#   t2 <- z0*wsb + (a[1]/(s[1]+1))*wsb^(s[1]+1)
#   t3 <- (z0+a[1]*wsb^s[1])*(wbf-wsb)
#   t4 <- (a[2]/(s[2]+1))*(wbf^(s[2]+1))
#   A <- t1 - t2 -(t3+t4)
#   return(A)
# }

nlsb_area3 <- function(wbf, WSEbf, z0, a, s, wsb) # pretty sure this is the right formula. Should run a test case or two.
{
  t1 <- WSEbf*wbf
  t2 <- z0*wsb + (a[1]/(s[1]+1))*wsb^(s[1]+1)
  
  t3_1 <- (z0+a[1]*wsb^s[1])*(wbf-wsb)
  t3_2 <- (a[2]/(s[2]+1))*(s[2]*wsb^(s[2]+1) + wbf*(wbf^s[2] - (s[2]+1)*wsb^s[2]))
  
  A <- t1 - t2 - (t3_1 + t3_2)
  return(A)
}

# Move toward estimating A0 in particular, since the A fluctuations are observed by SWOT.


