fit_shape <- function(dat)
{
  # Fitting cross sections to transects extracted using auto_transects.R
  # This algorithm works pretty well, especially if the cross-sections have "nice" shapes.
  # Use if multiplicative error is assumed
  
  # Fit an interpolating function to the smoothed data
  f1 <- approxfun(dat$width, y = dat$depth, method="linear",
                  yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)
  
  # Enter in with a depth and retrieve a width:
  d0 <- min(dat$depth)
  dn <- max(dat$depth) # bankfull, should correspond to bankfull width
  d.in <- seq(d0, dn, length.out = 100)
  nn<- length(d.in)
  
  # Need to revise to handle more cases
  flow_width <- vector(length = nn)
  rts.debug <- vector("list", length = nn)
  for (k in 1:nn)
  {
    realistic_wd <- function(x, d = d.in[k])
    {
      pred <- f1(x)-d
      return(pred)
    }
    rts <- tryCatch(uniroot.all(realistic_wd, c(0,1)), error = function(e) {NA})
    rts.debug[[k]] <- rts
    
    # Handle different cases (0, 1, 2, or more roots)
    if (length(rts)==0)
    {
      flow_width[k] <- 1
    }else if (length(rts)==1)
    {
      if (rts>0.5) {flow_width[k] <- 1-rts}
      else {flow_width[k] <- rts}
    }else if (length(rts)==2)
    {
      flow_width[k] <- rts[2] - rts[1]
    }else if (length(rts)%%2 == 0) # even, more than 2
    {
      sum1 <- 0
      N <- length(rts)
      for (i in seq(1, N, by = 2))
      {
        sum1 <- sum1 + rts[i+1] - rts[i]
      }
      flow_width[k] <- sum1
    }else # odd, more than 2
    {
      flow_width[k] <- NA # should handle this better
    }
    
  }
  
  # Set known points
  flow_width[1] <- 0
  flow_width[nn] <- 1
  
  # Fit the flow-width/depth relationship
  
  datfd <- data.frame(w = flow_width, d = d.in)
  
  datfd.no.zero <- datfd[-1,]
  datfd.no.zero <- na.omit(datfd.no.zero)
  
  fit.log <- lm(log(w)~log(d)-1, datfd.no.zero)
  
  s <- 1/fit.log$coefficients
  names(s) <- "shape"
  
  ##########
  plot.flag = TRUE
  # Plot, optionally
  if (plot.flag == TRUE)
  {
    x1 <- 0
    x2 <- 1
    xx <- seq(x1,x2,length.out = 100)
    
    par(mfrow=c(1,2), mar = c(5,5,5,5))
    plot(xx, f1(xx), type="l", xlim = c(0,1), ylim = c(min(dat$depth),1), 
         main = "", xlab = "normalized width", ylab = "normalized depth")
    points(dat)
    legend("top", legend = c("observed","fitted"), lty = c(NA, 1), pch = c(1, NA))
    
    plot(d.in, flow_width, xlab = "normalized depth", ylab = "normalized width",
         main = paste("Shape parameter = ", round(s,2)), cex = 1, col = "cyan")
    lines(datfd.no.zero$d, exp(predict(fit.log)))
  }
  ##########
  
  return(s)
}
  
# -------------------------------------------------------------------------------------------

fit_shape2 <- function(dat)
{
  # Revised version of fit_shape that doesn't muck around with the data as much
  # Created June 28, 2018. In development.
  
  # Fit an interpolating function to the smoothed data
  f1 <- approxfun(dat$width, y = dat$depth, method="linear",
                  yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)
  
  # Enter in with a depth and retrieve a width:
  d0 <- min(dat$depth)
  dn <- max(dat$depth) # bankfull, should correspond to bankfull width
  d.in <- seq(d0, dn, length.out = 100)
  nn<- length(d.in)
  
  # Need to revise to handle more cases
  flow_width <- vector(length = nn)
  rts.debug <- vector("list", length = nn)
  for (k in 1:nn)
  {
    realistic_wd <- function(x, d = d.in[k])
    {
      pred <- f1(x)-d
      return(pred)
    }
    rts <- tryCatch(uniroot.all(realistic_wd, c(0,1)), error = function(e) {NA})
    rts.debug[[k]] <- rts
    
    # Handle different cases (0, 1, 2, or more roots)
    if (length(rts)==0)
    {
      flow_width[k] <- 1
    }else if (length(rts)==1)
    {
      if (rts>0.5) {flow_width[k] <- 1-rts}
      else {flow_width[k] <- rts}
    }else if (length(rts)==2)
    {
      flow_width[k] <- rts[2] - rts[1]
    }else if (length(rts)%%2 == 0) # even, more than 2
    {
      sum1 <- 0
      N <- length(rts)
      for (i in seq(1, N, by = 2))
      {
        sum1 <- sum1 + rts[i+1] - rts[i]
      }
      flow_width[k] <- sum1
    }else # odd, more than 2
    {
      flow_width[k] <- NA # should handle this better
    }
    
  }
  
  # Fit the flow-width/depth relationship
  
  datfd <- data.frame(w = flow_width, d = d.in)
  datfd <- datfd[c(-1,-nn),] # remove the first and last points
  fit.log <- lm(log(w)~log(d)-1, datfd)
  
  s <- 1/fit.log$coefficients
  names(s) <- "shape"
  
  # Plot, optionally
  plot.flag = TRUE
  if (plot.flag == TRUE)
  {
    x1 <- 0
    x2 <- 1
    xx <- seq(x1,x2,length.out = 100)
    
    par(mfrow=c(1,2), mar = c(5,5,5,5))
    plot(xx, f1(xx), type="l", xlim = c(0,1), ylim = c(min(dat$depth),1), 
         main = "", xlab = "normalized width", ylab = "normalized depth")
    points(dat)
    legend("top", legend = c("observed","fitted"), lty = c(NA, 1), pch = c(1, NA))
    
    plot(w~d, data = datfd, xlab = "normalized depth", ylab = "normalized width",
         main = paste("Shape parameter = ", round(s,2)), cex = 1, col = "cyan")
    lines(exp(predict(fit.log)), datfd$d)
    lines(datfd$d, exp(predict(fit.log)))
  }
  
  return(s)
}

# -------------------------------------------------------------------------------------------
# Fitting by trial and error

#plot(w~d, datfd)
#dd <- seq(0,1,length.out = 30)
#s <- 1.5
#ww <- dd^(1/s)
#lines(dd,ww, type="l")

# -------------------------------------------------------------------------------------------

dingmansfunctions <- function()
{
  # Trying out some methods from Dingman's textbook - not really working right now - June 29, 2018
  
  # Dingman method 2
  
  y.bf <- 2.45
  psi.bf <- 5.77
  r.hat <- (y.bf/psi.bf)/(1-(y.bf/psi.bf))
  print(r.hat)
  
  # Dingman method 3a
  # Assumes symmetric cross section
  ll <- floor(length(dat$depth)/2)
  datf1 <- data.frame(bed = dat$depth, x = dat$width)
  datf1 <- datf1[1:ll,]
  # Remove zero values
  omit.ind <- which(datf1$bed==0)
  datf1 <- datf1[-omit.ind,]
  plot(bed~x, datf1)
  fit1 <- lm(log(bed)~log(x), datf1)
  summary(fit1)
  lines(datf1$x, exp(predict(fit1)))
  r <- fit1$coefficients[2]
  
  # Dingman method 3b
  r <- (fit1$coefficients[1]-log(bbf[[seg]]))/(log(2/wbf[[seg]]))
}

# ------------------------------------------------------------------------------------------

fit_nls <- function(dat)
{
  # nonlinear least squares
  # Use if additive error is assumed
  
  # Fit an interpolating function to the smoothed data
  f1 <- approxfun(dat$width, y = dat$depth, method="linear",
                  yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)
  
  # Enter in with a depth and retrieve a width:
  d0 <- min(dat$depth)
  dn <- max(dat$depth) # bankfull, should correspond to bankfull width
  d.in <- seq(d0, dn, length.out = 100)
  nn<- length(d.in)
  
  # Need to revise to handle more cases
  flow_width <- vector(length = nn)
  rts.debug <- vector("list", length = nn)
  for (k in 1:nn)
  {
    realistic_wd <- function(x, d = d.in[k])
    {
      pred <- f1(x)-d
      return(pred)
    }
    rts <- tryCatch(uniroot.all(realistic_wd, c(0,1)), error = function(e) {NA})
    rts.debug[[k]] <- rts
    
    # Handle different cases (0, 1, 2, or more roots)
    if (length(rts)==0)
    {
      flow_width[k] <- 1
    }else if (length(rts)==1)
    {
      if (rts>0.5) {flow_width[k] <- 1-rts}
      else {flow_width[k] <- rts}
    }else if (length(rts)==2)
    {
      flow_width[k] <- rts[2] - rts[1]
    }else if (length(rts)%%2 == 0) # even, more than 2
    {
      sum1 <- 0
      N <- length(rts)
      for (i in seq(1, N, by = 2))
      {
        sum1 <- sum1 + rts[i+1] - rts[i]
      }
      flow_width[k] <- sum1
    }else # odd, more than 2
    {
      flow_width[k] <- NA # should handle this better
    }
    
  }
  
  # Set known points
  flow_width[1] <- 0
  flow_width[nn] <- 1
  
  # Fit the flow-width/depth relationship
  
  datfd <- data.frame(w = flow_width, d = d.in)
  fit.nls <- nls(w~d^(1/s),start = list(s = 1),data=datfd)
  s <- coef(fit.nls)
  
  ##########
  plot.flag = TRUE
  # Plot, optionally
  if (plot.flag == TRUE)
  {
    x1 <- 0
    x2 <- 1
    xx <- seq(x1,x2,length.out = 100)
    
    par(mfrow=c(1,2), mar = c(5,5,5,5))
    plot(xx, f1(xx), type="l", xlim = c(0,1), ylim = c(min(dat$depth),1), 
         xlab = "normalized width", ylab = "normalized depth")
    points(dat)
    legend("top", legend = c("observed","fitted"), lty = c(NA, 1), pch = c(1, NA))
    
    plot(d.in, flow_width, xlab = "normalized depth", ylab = "normalized width",
         main = paste("Shape parameter = ", round(s,2)), cex = 1, col = "cyan")
    lines(datfd$d, predict(fit.nls))
  }
  ##########
  
  return(s)
  
}




