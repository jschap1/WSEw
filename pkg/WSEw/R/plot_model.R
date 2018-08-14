#' Plot modeled cross section
#' 
#' @param model
#' @param type
#' @param WSEw required for nonlinear models
#' @export

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
    warning("sbm not yet implemented")
    
    
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
  
  if (add)
  {
    lines(x,b, ...)
  } else
  {
    plot(x,b, type = "l", ...)
  }
  
}

#' Calculate b(x) for the slope break model
#' @export
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


#' Calculate b(x) for the nonlinear slope break model
#' @export
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
