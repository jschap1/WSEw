fit_nls2 <- function(dat)
{
  # nonlinear least squares
  # Use if additive error is assumed
  # Uses two functions to fit the data
  
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
  
  # Split the dataset at the slope break
  plot(flow_width)
  plot(diff(flow_width))
  plot(diff(flow_width,differences = 2), ylim = c(-0.01, 0.01))
  abline(0,0)
  
  # Fit the flow-width/depth relationship
  
  datfd <- data.frame(w = flow_width, d = d.in)
  fit.nls <- nls(w~d^(1/s),start = list(s = 1),data=datfd)
  s <- coef(fit.nls)
  
  # Test for structural change
  testres <- sctest(fit.nls)
  testres <- sctest(fit.nls, functional = "supLM")
  testres # returns f(efp) and p-value
  # null hypothesis is parameter stability
  # alternative is that the value of the parameter changes
  # reject the null hypothesis if the p-value is small
  # I don't think this is working so well.
  # Stick to Mersel's "easy" method
  # Will need to specify shape parameter a priori.
  
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