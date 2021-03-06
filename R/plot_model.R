#' Plot modeled cross section
#' 
#' Computes modeled bed elevation and plots the cross section
#' @export
#' @param model
#' @param type Either "linear", "sb", "sbm", "nl", or "nlsb"
#' @param WSEw required for nonlinear models
#' @param add TRUE to add plot to current figure, FALSE to create a new figure. Default is FALSE
#' @param ... additional parameters for plot or lines
#' @details Eventually, it might make sense for this to return the x and b values in list similar to cross_sections from auto_transects
#' @examples 
#' plot(cross_sections$x[[2]], cross_sections$b[[2]], type = "l")
#' plot_model(lf1, type = "linear", col = "red", add = TRUE)

plot_model <- function(model, type, WSEw = NULL, add = FALSE, ...)
{
  if (type == "linear")
  {
    
    x0 <- 0
    xf <- max(model$model$w)
    z0 <- as.numeric(coef(model)[1])
    s <- as.numeric(coef(model)[2])
    x <- seq(x0, xf, length.out = 100)
    b <- z0 + abs(2*x - (x0+xf))*s
    
    } else if (type == "sb")
  {
    
    m1 <- model[[1]]
    m2 <- model[[2]]

    s1 <- as.numeric(coef(m1)[2])
    s2 <- as.numeric(coef(m2)[1])
    wsb <- max(m1$model$w)
    wbf <- max(m2$model$`I(w - wb)`) + wsb # accounting for the offset
    x0 <- 0
    xf <- wbf
    x <- seq(x0, xf, length.out = 100)
    z0 <- predict(m1, newdata = data.frame(w=0))
    b <- sb_model_xs(x, z0, s1, s2, wsb, x0, xf)
    
  } else if (type == "sbm")
  {
    warning("sbm only implemented for the case of exactly 3 slope breaks")
    
    m <- model[[1]]
    wbf <- max(m$model$w)
    x0 <- 0
    xf <- wbf
    
    cc <- as.numeric(coef(m))
    cc <- cc[which(cc!=0)]
    alpha1 <- cc[1] # intercept term; see derivation for notation here
    u <- cc[-1]
    beta <- 2*cumsum(u)
    wsb <- as.numeric(m$psi[,1])
    xsb <- wsb/2
    
    x01 <- seq(x0, xsb[1], length.out = 30)
    x12 <- seq(xsb[1], xsb[2], length.out = 30)
    x23 <- seq(xsb[2], xsb[3], length.out = 30)
    x34 <- seq(xsb[3], xf/2, length.out = 30)
    
    b01 <- alpha1 + beta[1]*x01
    b12 <- alpha1 + beta[1]*xsb[1] + beta[2]*(x12-xsb[1])
    b23 <- alpha1 + beta[1]*xsb[1] + beta[2]*(xsb[2]-xsb[1]) + beta[3]*(x23-xsb[2])
    b34 <- alpha1 + beta[1]*xsb[1] + beta[2]*(xsb[2]-xsb[1]) + beta[3]*(xsb[3] - xsb[2]) + beta[4]*(x34-xsb[3])
    
    xr <- c(x01, x12, x23, x34) # right
    b <- c(b01, b12, b23, b34)

    xl <- -xr # left
    xr <- xr + (xf-x0)/2
    xl <- xl + (xf-x0)/2
  
  } else if (type == "nl")
  {
    
    x0 <- 0
    xf <- max(WSEw$w)
    z0 <- as.numeric(coef(model)[1])
    a <- as.numeric(coef(model)[2])
    s <- as.numeric(coef(model)[3])
    x <- seq(x0, xf, length.out = 100)
    b <- z0 + a*abs(2*x-x0-xf)^s
      
  } else if (type == "nlsb")
  {

    m1 <- model[[1]]
    m2 <- model[[2]]
    
    z0 <- predict(m1, newdata = data.frame(w=0))
    a1 <- as.numeric(coef(m1)[2])
    a2 <- as.numeric(coef(m2)[1])
    s1 <- as.numeric(coef(m1)[3])
    s2 <- as.numeric(coef(m2)[2])
    
    sb.ind <- as.numeric(attributes(model)) # index of slope break
    wsb <- WSEw$w[sb.ind] # width at slope break
    wbf <- max(WSEw$w)
    
    x0 <- 0
    xf <- wbf
    x <- seq(x0, xf, length.out = 100)
    b <- nlsb_model_xs(x, z0, a1, a2, s1, s2, wsb, x0, xf)
  }
  
  if (type == "sbm")
  {
    if (add)
    {
      lines(xl, b, type = "l", ...)
      lines(xr, b)
    } else
    {
      plot(xl, b, type = "l", ...)
      lines(xr, b)
    }
    
  }
  else
  {
    if (add)
    {
      lines(x,b, ...)
    } else
    {
      plot(x,b, type = "l", ...)
    }
  }

}

# --------------------------------------------------------------------------------------------------------------------


#' Calculate slope break model cross section
#' 
#' Calculate b(x) for the slope break model
#' @export
#' @param x vector of x coordinates
#' @param z0 minimum bed elevation
#' @param s1 slope parameter of the first linear model
#' @param s2 slope parameter of the second linear model
#' @param wsb width at the breakpoint
#' @param x0 minimum x value, default is 0
#' @param xf maximum x value
#' @details Runs in a loop, could be vectorized for speed
#' @return b modeled bed elevation
#' @example b <- sb_model_xs(x, z0, s1, s2, wsb, x0, xf)

sb_model_xs <- function(x, z0, s1, s2, wsb, x0, xf)
{
  # x coordinates of the beginning and end of the lower part of the cross section
  xsb1 <- 0.5*(x0+xf-wsb) 
  xsb2 <- 0.5*(x0+xf+wsb)
  
  n <- length(x)
  b <- vector(length = n)
  for (i in 1:n)
  {
    if (x[i] <= xsb1)
    {
      b[i] <- z0 + s1*(wsb) + 2*s2*(xsb1 - x[i])
    } else if (x[i]> xsb1 & x[i]<xsb2)
    {
      b[i] <- z0 + s1*abs(xsb2+xsb1-2*x[i])
    } else if (x[i]>xsb2)
    {
      b[i] <- z0 + s1*(wsb) + 2*s2*(x[i]-xsb2)
    }else
    {
      b[i] <- NA
    }
  }
  return(b)
}



# --------------------------------------------------------------------------------------------------------------------

#' Calculate NLSB model cross section
#'
#' Calculate b(x) for the nonlinear slope break model
#' @export
#' @param x vector of x coordinates
#' @param z0 minimum bed elevation
#' @param a1 multiplicative parameter of the first nonlinear model
#' @param a2 multiplicative parameter of the second nonlinear model
#' @param s1 exponent of the first nonlinear model
#' @param s2 exponent of the second nonlinear model
#' @param wsb width at the breakpoint
#' @param x0 minimum x value, default is 0
#' @param xf maximum x value
#' @details Runs in a loop, could be vectorized for speed
#' @return b modeled bed elevation
#' @example b <- nlsb_model_xs(x, z0, a1, a2, s1, s2, wsb, x0, xf)

nlsb_model_xs <- function(x, z0, a1, a2, s1, s2, wsb, x0, xf)
{
  # x coordinates of the beginning and end of the lower part of the cross section
  xsb1 <- 0.5*(x0+xf-wsb) 
  xsb2 <- 0.5*(x0+xf+wsb)
  
  n <- length(x)
  b <- vector(length = n)
  for (i in 1:n)
  {
    if (x[i]> xsb1 & x[i]<xsb2) 
    {
      # below slope break
      b[i] <- z0 + a1*abs(x0+xf-2*x[i])^s1
    }else 
    {
      # above slope break
      b[i] <- z0 + a1*wsb^s1 + a2*(abs(x0+xf-2*x[i])^s2 - wsb^s2)
    }
  }
  return(b)
}
