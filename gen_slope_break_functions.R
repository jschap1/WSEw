# Generalized slope break method supporting functions
# July 17, 2018


get_width2 <- function(WSE, f1, x.max, delx = 1)
{
  # Numerically estimates the top width
  
  n <- length(WSE)
  flow_width <- vector(length = n)
  
  b.trap <- f1(seq(0, x.max, by = delx))
  
  for (k in 1:n)
  {
    flow_width[k] <- delx*sum(WSE[k]>b.trap)
    
    # Once the WSE is higher than any part of the bed elevation, stop
    if (sum(WSE[k] > max(b.trap)))
    {
      flow_width[k] <- NA
    }
  }
  
  return(flow_width)
}

# --------------------------------------------------------------------

fit_nls2 <- function(dat)
{
  # Better version than fit_nls
  
  
  
  return(s)
}

# --------------------------------------------------------------------

get_width <- function(WSE, f1, x.max)
  # Use get_width2, it is more robust
{
  n <- length(WSE)
  flow_width <- vector(length = n)
  rts.debug <- vector("list", length = n)
  
  f1.min <- min(f1(seq(0, x.max)))
  f1.max <- max(f1(seq(0, x.max)))
  
  for (k in 1:n)
  {
    realistic_wd <- function(x, d = WSE[k])
    {
      pred <- f1(x)-d
      return(pred)
    }
    rts <- tryCatch(uniroot.all(realistic_wd, c(0, x.max)), error = function(e) {NA})
    rts.debug[[k]] <- rts
    
    # Handle different cases (0, 1, 2, or more roots)
    if (length(rts)==0)
    {
      if (WSE < f1.min)
      {
        flow_width[k] <- 0
      } else if (WSE[k] > f1.max)
      {
        flow_width[k] <- x.max
      } else
      {
        flow_width[k] <- NA
      }
    }
    
    else if (length(rts)==1) # figure this one out...
    {
      if (rts>x.max/2) 
      {
        if (WSE[k] > f1.max)
        {
          flow_width[k] <- x.max-rts
        } else 
        {
          flow_width[k] <- 0
        }

      }
      else
      {
        flow_width[k] <- rts
      }
    } 
    
    else if (length(rts)==2)
    {
      # estimate derivative at left bank
      df1 <- (f1(rts[1]+1) - f1(rts[1]-1))/2
      if (df1 > 0)
      {
        rts <- c(rts[-1], x.max)
      } 
      
      flow_width[k] <- rts[2] - rts[1]
    }
    
    else if (length(rts)%%2 == 0) # even, more than 2
    {
      # estimate derivative at left bank
      df1 <- (f1(rts[1]+1) - f1(rts[1]-1))/2
      if (df1 > 0)
      {
        rts <- c(rts[-1], x.max)
      }
      
      sum1 <- 0
      N <- length(rts)
      for (i in seq(1, N, by = 2))
      {
        sum1 <- sum1 + rts[i+1] - rts[i]
      }
      flow_width[k] <- sum1
    }
  
    else # odd, more than 2
    {
      # estimate derivative at left bank
      df1 <- (f1(rts[1]+1) - f1(rts[1]-1))/2
      if (df1 > 0)
      {
        rts <- c(0, rts)
      } else if (df1 <= 0)
      {
        rts <- c(rts, x.max)
      }
      
      sum1 <- 0
      N <- length(rts)
      for (i in seq(1, N, by = 2))
      {
        sum1 <- sum1 + rts[i+1] - rts[i]
      }
      flow_width[k] <- sum1
    }
    
  }
  
  # Set known points
  #flow_width[1] <- 0
  #flow_width[n] <- x.max
  
  return(flow_width)
}




# ----------------------------------------------------------------------------------------------------




calc_width <- function(dat)
{
  # Calculates the width for a given depth
  
  # Fit an interpolating function to the channel geometry
  f1 <- approxfun(dat$x, y = dat$d, method="linear",
                  yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)
  
  
}

get_w_WSE_relation <- function(dat)
{  
  # Returns datfd, flow width vs. flow depth
  
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
  
  datfd <- data.frame(w = flow_width, d = d.in)
  return(datfd)
  
}


# ------------------------------------------------------------------------------------------------

test_method <- function(datfd)
{
  # Testing for linear vs. slope break method
  nn <- dim(datfd)[1]
  
  # Estimate first derivative using forward difference
  fd <- vector(length = nn)
  for (k in 1:(nn-1))
  {
    fd[k] <- (datfd$d[k+1] - datfd$d[k])/(datfd$w[k+1] - datfd$w[k])
  }
  
  thres <- 0.015 # this is way too small, given that I've normalized the widths and depths. 
  if (max(fd) - min(fd) < thres)
  {
    t <- 0 # "linear"
  }else
  {
    t <- 1 # "slope break"
  }
  
  # (datfd$d[2] - datfd$d[1])/(datfd$w[2] - datfd$w[1])
  
  # mean(fd) # mean slope of the width-depth relationship
  # compute WSE_min by extrapolation...
  
  return(t) # 0 for linear, 1 for slope break
}

}




# ------------------------------------------------------------------------------------------------
# Scrap

  # nonlinear least squares
  # Use if additive error is assumed
  # Uses two functions to fit the data
  

  
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
  





