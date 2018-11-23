#' Calculate modeled bankfull flow area
#' 
#' Assumes some model, such as linear, nonlinear, or slope break for the channel cross section
#' @export
#' @param model fitted model to the WSE-w data
#' @param type type of model. Default is linear. Options are linear, sb, sbm, nl, and nlsb
#' @param WSEw WSE-w data for the fit. Only required if fit is nonlinear
#' @return A cross-sectional area of flow
#' @examples 
#' A.l <- calc_modelA(lf[[1]], type = "linear")
#' A.sb <- calc_modelA(sbm[[1]], type = "sbm")

calc_model_A <- function(model, type, WSEw = NULL)
{
  
  if (type == "linear")
  {
    
    z0 <- coef(model)[1]
    a <- coef(model)[2]
    h_obs <- model$model$WSE
    hbf <- max(h_obs)
    wbf <- (hbf - z0)/a
    A <- 0.5*wbf*(hbf - z0)
    
  }else if (type == "sb") # it would be more straightforward to use segmented for everything, including the single SB model
  {
    m1 <- model[[1]] # lower model
    m2 <- model[[2]] # upper model 
    wsb <- max(m1$model$w)
    wbf <- max(m2$model$`I(w - wb)`) + wsb # accounting for the offset
    WSEsb <- max(fitted(m1))
    hbf <- max(fitted(m2))
    z0 <- as.numeric(coef(m1)[1])
    A <- sb_area(wsb, WSEsb, z0, hbf, wbf)   
    
  }else if (type == "sbm")
  {
    
    model <- model[[1]]
    
    wbf <- max(model$model$w)
    # hbf <- max(model$model$WSE)
    hbf <- max(fitted(model))
    
    w.brk <- as.numeric(summary(model)$psi[,2]) # width breakpoints, in increasing order
    WSE.brk <- predict(model, newdata = data.frame(w = w.brk))
    
    w <- c(w.brk, wbf)
    WSE <- c(WSE.brk, hbf)
    # z0 <- min(model$model$WSE)
    # z0 <- predict(model, newdata = data.frame(w = 0))
    z0 <- as.numeric(coef(model)[1])
    
    A <- sbm_area(w, WSE, z0) 
    
  }else if (type == "nl")
  {
    
    z0 <- model$z0
    a <- model$a
    s <- model$s
    h_obs <- WSEw$WSE
    hbf <- max(h_obs)
    wbf <- ((hbf - z0)/a)^(1/s)
    A <- nl_area(wbf, hbf, z0, a, s)
    
  }else if (type == "nlsb")
  {

    a1 <- model$a1
    a2 <- model$a2
    s1 <- model$s1
    s2 <- model$s2
    z0 <- model$z0

    h_obs <- WSEw$WSE
    hbf <- max(h_obs)
    wbf <- ((hbf - z0)/a1)^(1/s1)
    
    sb.ind <- model$sb.ind # index of slope break
    wb <- WSEw$w[sb.ind] # width at slope break

    # if a2 is negative, ignore the cross sectional area above the slope break;
    # it cannot be estimated well, given the measurement error
    # if (a2 < 0) {a2 <- 0}
    
    A <- nlsb_area4(wbf, hbf, z0, a1, a2, s1, s2, wb)
    
    # numerical method (alternative, doesn't work perfectly...)
    # w.1 <- seq(0, wb, length = 1e3) # arbitrary (but large) length
    # w.2 <- seq(wb, wbf, length = 1e3)
    # h.1 <- z0 + a1*w.1^s1
    # h.2 <- z0 + a1*wb^s1 + a2*(w.2^s2 - wb^s2) 
    # WSEw1 <- data.frame(WSE = c(h.1, h.2), w = c(w.1, w.2))
    # A <- calc_A_from_WSEw(WSEw1)
    
  }
  
  return(A)
}

# r=7, k=5, m = 1
# A.nlsb[7,5,1]

# --------------------------------------------------------------------------------------------------

#' Calculates area of a slope break cross section
#' @export
 
sb_area <- function(wsb, WSEsb, z0, hbf, wbf)   
{
  A <- 0.5*(WSEsb-z0)*wsb + 0.5*(wsb+wbf)*(hbf - WSEsb)
  return(A)
}

# --------------------------------------------------------------------------------------------------

#' Calculates area of a multi-slope break cross section
#' @param w width
#' @param WSE water surface elevation. The lowest value is z0 and the highest is bankfull WSE.
#' @export

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
#' @export

nl_area <- function(wbf, hbf, z0, a, s)
{
  A <- (hbf - z0)*wbf - (a/(s+1))*wbf^(s+1)
  # A <- (hbf - z0)*wbf - (a/(s+1))*wbf^(s+1)
  
  return(A)
}

# --------------------------------------------------------------------------------------------------

#' Calculate area of a nonlinear slope break cross section
#' @param s shape parameters for nonlinear fits (the exponents)
#' @param a multiplicative parameters for the nonlinear fits
#' @export

# Move toward estimating A0 in particular, since the A fluctuations are observed by SWOT.

nlsb_area4 <- function(wbf, hbf, z0, a1, a2, s1, s2, wb)
{
  
  t1 <- (hbf - z0)*wb
  t2 <- (a1/(s1+1))*wb^(s1+1)
  t3 <- (hbf - z0 - a1*wb^s1)*(wbf-wb)
  t4 <- (a2*wb^s2)*(wbf-wb) # this term can get very large if a2 is not small
  t5 <- (a2/(s2+1))*(wbf-wb)^(s2+1)
  
  A <- t1 - t2 + t3 + t4 - t5
  
  return(A)
}

