#' Calculate modeled wetted perimeter
#' 
#' Calculates wetted perimeter for a modeled cross section using one of several types of geometric channel models
#' @export
#' @param model model object
#' @param w flow width vector, needed for nonlinear models
#' @details Uses correct formulas to calculate WP for the modeled cross sections. 
#' For all the UMRB cross sections, width is much larger than depth, so WP = w, though.
#' @return WP wetted perimeter
#' @example WP <- calc_model_WP(lf1, type = "linear")

calc_model_WP <- function(model, type, w = NULL)
{
   if (type == "linear")
   {
   
     wbf <- max(model$model$w)
	   z0 <- predict(model, newdata = data.frame(w=0))
	   WSEbf <- max(fitted(model))
     WP <- calc_WP_linear(wbf, WSEbf, z0)
   
   } else if (type == "sb")
   {
   
     m1 <- model[[1]]
	   m2 <- model[[2]]
     wsb <- max(m1$model$w)
     wbf <- max(m2$model$`I(w - wb)`) + wsb
	   z0 <- predict(m1, newdata = data.frame(w=0))
     WSEsb <- max(fitted(m1))
	   WSEbf <- max(fitted(m2))
     WP <- calc_WP_sb(wbf, WSEbf, z0, wsb, WSEsb)
   
   } else if (type == "sbm")
   {
   
     model <- model[[1]] # needs work, see area calculations code

     wbf <- max(model$model$w)
     WSEbf <- max(fitted(model))
     
     w.brk <- as.numeric(summary(model)$psi[,2]) # width breakpoints, in increasing order
     WSE.brk <- predict(model, newdata = data.frame(w = w.brk))
     
     z0 <- predict(model, newdata = data.frame(w = 0))
     
     w <- c(0, w.brk, wbf)
     WSE <- c(z0, WSE.brk, WSEbf)

     WP <- calc_WP_sbm(w, WSE)

   } else if (type == "nl")
   {
     
     wbf <- max(w)
     a <- as.numeric(coef(model)[2])
     s <- as.numeric(coef(model)[3])
     WP <- calc_WP_nl(wbf, a, s)
   
   } else if (type == "nlsb")
   {
   
     m1 <- model[[1]] # lower model
     m2 <- model[[2]] # upper model 
     
     wbf <- max(w)
     sb.ind <- as.numeric(attributes(model)) # index of slope break
     wsb <- w[sb.ind] # width at slope break
     
     a1 <- as.numeric(coef(m1)[2])
     a2 <- as.numeric(coef(m2)[1])
     s1 <- as.numeric(coef(m1)[3])
     s2 <- as.numeric(coef(m2)[2])

     WP <- calc_WP_nlsb(wbf, wsb, a1, a2, s1, s2)
   
   }
  return(WP)
}

# ---------------------------------------------------------------------------

#' @export
calc_WP_linear <- function(wbf, WSEbf, z0)
{
  WP <- (wbf^2 + 4*(WSEbf - z0)^2)^0.5
  return(WP)
}

# ---------------------------------------------------------------------------

#' @export
calc_WP_sb <- function(wbf, WSEbf, z0, wsb, WSEsb)
{
  WP <- (wsb^2 + 4*(WSEsb - z0)^2)^0.5 + ((wbf-wsb)^2 + 4*(WSEbf - WSEsb)^2)^0.5
  return(WP)
}

# ---------------------------------------------------------------------------

calc_WP_sbm <- function(w, WSE)
{
  WP <- 0
  n <- length(w)
  for (k in 1:(n-1))
  {
    WP <- WP + ((w[k+1]-w[k])^2 + 4*(WSE[k+1] - WSE[k])^2)^0.5
  }
  return(WP)
}

# ---------------------------------------------------------------------------

calc_WP_nl <- function(wbf, a, s)
{
  # uses a built in numerical integration function in R
  # there is also a pracma alternative "integral" with more options
  fun <- function(w, a, s)
  {
    f <- (1+(2*a*(s-1)*w^(s-1))^2)^0.5
    return(f)
  }
  # print(paste("f", fun))
  # print(paste("wbf", wbf))
  # print(paste("a", a))
  # print(paste("s", s))
  if (s<=0.1) # does not work if s -> 0 
  {
    WP <- NA
    return(WP)
  }
  WP <- integrate(f = fun, lower = 0, upper = wbf, a = a, s = s) # passing additional arguments to fun 
  WP <- WP$value
  return(WP)
}

# ---------------------------------------------------------------------------

calc_WP_nlsb <- function(wbf, wsb, a1, a2, s1, s2)
{

  if (s1 <= 0 | s2 <= 0) # does not work if s<=0 
  {
    WP <- NA
    return(WP)
  }
  
  fun1 <- function(w, a1, s1)
  {
    f <- (1+(2*a1*(s1-1)*w^(s1-1))^2)^0.5
    return(f)
  }
  
  fun2 <- function(w, a2, s2)
  {
    f <- (1+(2*a2*(s2-1)*w^(s2-1))^2)^0.5
    return(f)
  }
  
  if (s1<=0.1 | s2 <= 0.1) # does not work if s -> 0 
  {
    WP <- NA
    return(WP)
  }
  
  WP1 <- integrate(f = fun1, lower = 0, upper = wsb, a1 = a1, s1 = s1)
  WP2 <- integrate(f = fun2, lower = wsb, upper = wbf, a2 = a2, s2 = s2)
  WP <- WP1$value + WP2$value
  
  return(WP)
}

