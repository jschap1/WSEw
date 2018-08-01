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
    
    # Once the WSE is higher than any other part of the bed elevation,
    # the width stops changing, as rectangular walls as assumed
    if (sum(WSE[k] > max(b.trap)))
    {
      flow_width[k] <- flow_width[k-1]
    }
  }
  
  return(flow_width)
}

# ------------------------------------------------------------------------------------------------

test_linear <- function(WSE, w, thres = 0.015)
{
  # Tests if the linear method is appropriate for depth estimation
  
  nn <- length(WSE)
  s1 <- vector(length = nn)
  for (k in 1:(nn-1))
  {
    s1[k] <- (WSE[k+1]-WSE[k])/(w[k+1]-w[k]) # forward difference
  }
  s1[nn] <- s1[nn-1]
  # remove infinite values (where width did not change)
  s1[s1==Inf] <- NA
  s1[s1==-Inf] <- NA
  
  maxdiff <- max(s1, na.rm = TRUE) - min(s1, na.rm=TRUE)
  
  if (maxdiff < thres)
  {
    print("Use the linear method")
    return(s1)
  } else
  {
    print("Do not use the linear method")
    return(NULL)
  }

}


# ------------------------------------------------------------------------------------------

test_slope_break <- function(WSE, w, thres = 0.015, window = 4, m = FALSE)
{
  # Tests whether the slope break method is appropriate for depth estimation
  # m = flag for multiple slope break method
  
  # calculate mean of the slope for the highest few water surfaces
  nn <- length(WSE)
  s1 <- vector(length = nn)
  for (k in 1:(nn-1))
  {
    s1[k] <- (WSE[k+1]-WSE[k])/(w[k+1]-w[k]) # forward difference
  }
  s1[nn] <- s1[nn-1]
  
  # deal with infinite values (where width did not change)
  s1[s1==Inf] <- 9999
  s1[s1==-Inf] <- -9999
  #s1[s1==Inf] <- NA
  #s1[s1==-Inf] <- NA
  
  s1.high <- s1[(nn-(window-1)):nn] 
  s.bar <- mean(s1.high, na.rm = TRUE)
  
  # compare to s1 at lower water surfaces to find a slope break
  for (j in 1:(nn-window-1))
  {
    #print(j) # debugging
    #print(length(s1))
    #print(nn-(window+j))
    #print(s1[nn-(window+j)])
    #print("s1")
    #print(s1[nn-(window+j)])
    #print("s.bar")
    #print(s.bar)
    if ((s1[nn-(window+j)] < 0.3*s.bar) | (s1[nn-(window+j)] > 1.7*s.bar))
    {
      sb.ind <- nn-(window+j)
      print(paste("There is a slope break at ", sb.ind))
      break
    } else
    {
      s.bar <- mean(s.bar, s1[nn-(window)-j], na.rm = TRUE)
      sb.ind <- FALSE
    }
  }
  
  if (sb.ind == FALSE)
  {
    print("Do not use the slope break method")
    return(NULL)
  }
  
  if (!m) # if using multi-slope break method, then report all the slope breaks
  {
    # Check that behavior is linear below the slope break
    s2 <- s1[1:sb.ind-1]
    maxdiff <- max(s2, na.rm = TRUE) - min(s1, na.rm=TRUE)
    if (maxdiff < thres)
    {
      print("Use the slope break method")
      return(sb.ind)
    } else
    {
      print("Do not use the slope break method")
      return(NULL)
    }
  } else
  {
    return(sb.ind)
  }
  
}

# ------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------

preprocess <- function(seg, w, wbf, d, dbf, rescale = FALSE)
{
  
  if (wbf[[seg]]==0)
  { # error handling
    return(NA)
  }
  w.fract <- w[[seg]]/wbf[seg]
  d.fract <- d[[seg]]/dbf[[seg]]
  
  dat <- data.frame(width=w.fract, depth=d.fract)
  
  if (rescale)
  { 
    # option to rescale the depths to 0-1
    num <- dat$depth - min(dat$depth, na.rm = T)
    denom <- max(dat$depth, na.rm = T) - min(dat$depth, na.rm = T)
    dat$depth <- num/denom
  }
  
  return(dat)
}

# ------------------------------------------------------------------------------------------

reach_avg <- function(xWSEw, l = 10000, res = 5)
{
  # Performs reach averaging on the WSE, w data
  #
  # INPUTS
  # l = reach length, m
  # res = cross-section spacing, m
  # WSEw pairs for each cross section: xWSEw
  #
  # OUTPUTS
  # Reach averaged WSE and w pairs: rWSEw
  
  pix <- l/res # the number of sections in a reach
  nseg <- length(xWSEw) # number of cross sections
  nn <- dim(xWSEw[[1]])[1]
  
  # Create matrices with WSE, w values
  WSEmat <- array(dim = c(nn, nseg))
  wmat <- array(dim = c(nn, nseg))
  for (seg in 1:nseg)
  {
    WSEmat[,seg] <- xWSEw[[seg]]$WSE
    wmat[,seg] <- xWSEw[[seg]]$w
  }
  
  # Do reach averaging
  nr <- nseg - (pix - 1) # number of reaches
  rWSEw <- vector(length = nr, "list")
  for (r in 1:nr)
  {
    rWSEw[[r]] <- data.frame(WSE = rowSums(WSEmat[,r:(r+(pix-1))], na.rm = TRUE)/pix
                             , w = rowSums(wmat[,r:(r+(pix-1))], na.rm = TRUE)/pix)
  }
  
  # Developing the averaging rule
  # rWSEw[[r]] <- data.frame(WSE = rowSums(WSEmat[,1:4])/pix # may want to handle NAs here
  #                         , w = rowSums(wmat[,1:4])/pix)

  # Note: the wbf values as used for defining channel exposure have changed. 
  # Now, it would make sense to use max(reach average w) as the wbf value.
  
  return(rWSEw)
}

# ------------------------------------------------------------------------------------------

xs_reach_avg <- function(cross.sections, l = 10000, res = 5)
{
  # Performs reach averaging on the observed bathymetry data
  # that is, directly on the cross-section geometry
  #
  # In development
  # Actually, it is not straightforward to compute a reach-average cross-section because 
  # the width of the cross section changes with distance downstream
  
  pix <- l/res # the number of sections in a reach
  nseg <- length(cross.sections) # number of cross sections

  # Create matrices with WSE, w values
  WSEmat <- array(dim = c(nn, nseg))
  wmat <- array(dim = c(nn, nseg))
  for (seg in 1:nseg)
  {
    xmat[,seg] <- cross.sections[[seg]]$x
    dmat[,seg] <- cross.sections[[seg]]$d
    bmat[,seg] <- cross.sections[[seg]]$b
  }
  
  # Do reach averaging
  nr <- nseg - (pix - 1) # number of reaches
  rCS <- vector(length = nr, "list") # reach averaged cross section geometry
  for (r in 1:nr)
  {
    rCS[[r]] <- data.frame(WSE = rowSums(WSEmat[,r:(r+(pix-1))], na.rm = TRUE)/pix
                             , w = rowSums(wmat[,r:(r+(pix-1))], na.rm = TRUE)/pix)
  }
  
  
}
  
# ------------------------------------------------------------------------------------------

calc_WSEw <- function(interval = 0.05, dx = 1)
{
  # Calculates WSE-w pairs for each cross section
  #
  # INPUTS
  # tdist, d, b, dbf, b.min for each segment
  # interval, number of points to sample in each cross-section
  # dx, used for computing width with trapezoidal rule
  #
  # OUTPUTS
  # WSE-w pairs for each cross section: xWSEw, cross section geometry
  
  ###########
  ########### loading/preprocessing that I would like to take out of the function...
  load("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Transects/pool_21_5m.rda")
  library(rootSolve)
  source("sb_functions.R")
  
  # Remove empty channels
  rm.ind <- which(channel.pix<=10)
  transects.depth <- transects.depth[-rm.ind]
  b <- b[-rm.ind]
  tdist <- tdist[-rm.ind]
  
  # Pre-processing
  nseg <- length(b)
  d <- transects.depth
  dbf <- unlist(lapply(d, max, na.rm=TRUE)) # bankfull depth (m)
  b.min <- unlist(lapply(b, min))
  ###########
  ###########
  
  cross.sections <- vector(length = nseg, "list")
  xWSEw <- vector(length = nseg, "list")
  for (seg in 1:nseg)
  {
    # Estimate WSE-w relationship
    cross.section <- data.frame(x = tdist[[seg]], d = d[[seg]], b = b[[seg]])
    cross.sections[[seg]] <- cross.section
    
    # Fit a linear spline to the observed channel geometry
    f1 <- approxfun(cross.section$x, y = cross.section$b, method="linear",
                    yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)
    
    # Make a vector of WSE(t) values to enter in with
    if (dbf[seg] < 0){next} # error check
    WSE <- seq(b.min[seg], b.min[seg] + dbf[seg], by = interval*dbf[seg])
    
    # Estimate the flow width for each WSE value
    w <- get_width2(WSE, f1, x.max = max(cross.section$x), delx = dx)
    
    # Put WSE-w info into a data frame
    WSEw <- data.frame(WSE = WSE, w = w)
    #WSEw <- na.omit(WSEw) # omit NAs, do not do this yet if reach averaging
    # really, should figure out how to do this without getting NA values
    
    xWSEw[[seg]] <- WSEw
  }
  
  # Remove null cross sections, those with no measurements (optional)
  rm.ind <- which(unlist(lapply(xWSEw, is.null)))
  xWSEw <- xWSEw[-rm.ind] # this doesn't always seem to work...
  
  result <- list(WSEw = xWSEw, xs = cross.sections)
  return(result)
}


# ------------------------------------------------------------------------------------------

check_below_sb <- function(WSE, w, sb.ind, thres = 0.015)
{
  # calculate mean of the slope for the highest few water surfaces
  nn <- length(WSE)
  s1 <- vector(length = nn)
  for (k in 1:(nn-1))
  {
    s1[k] <- (WSE[k+1]-WSE[k])/(w[k+1]-w[k]) # forward difference
  }
  s1[nn] <- s1[nn-1]
  
  # deal with infinite values (where width did not change)
  s1[s1==Inf] <- 9999
  s1[s1==-Inf] <- -9999
  
  # Check that behavior is linear below the slope break
  s2 <- s1[1:sb.ind-1]
  maxdiff <- max(s2, na.rm = TRUE) - min(s1, na.rm=TRUE)
  if (maxdiff < thres)
  {
    print("Use the slope break method")
    return(sb.ind)
  } else
  {
    print("Do not use the slope break method")
    return(NULL)
  }
}

# ------------------------------------------------------------------------------------------

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

test_gen_sb <- function(d.in, w, sb.thres = 0.3, window = 4)
{
  # Test for shape break aka generalized slope break
  # criterion: does s parameter change by a lot as you fit successive sets of points?
  
  nn <- length(d.in)
  high.ind <- (nn-(window-1)):nn
  datfd <- data.frame(w = w[high.ind], d = d.in[high.ind])
  fit.nls <- nls(w~d^(1/s),start = list(s = 1), data=datfd)
  s.high <- coef(fit.nls)
  sb.ind <- NULL

  # compare to s at lower water surfaces to find a "slope break"
  for (j in 1:(nn-2*window))
  {
    low.ind <- (nn-(2*window-1+j)):(nn-(window+j))
    datfd.low <- data.frame(w = w[low.ind], d = d.in[low.ind]) 
    fit.nls <- nls(w~d^(1/s),start = list(s = 1), data=datfd.low)
    s.low <- coef(fit.nls) 
    
    # Debugging
    #print(s.low)
    #print(s.high)
    
    if (s.low < sb.thres*s.high)
    {
      sb.ind <- low.ind[window]
      print(paste("There is a shape break at", low.ind[window]))
      break
    } 
    
    else
    {
      high.ind <- (nn-(window-1+j)):(nn)
      datfd <- data.frame(w = w[high.ind], d = d.in[high.ind])
      fit.nls <- nls(w~d^(1/s),start = list(s = 1), data=datfd)
      s.high <- coef(fit.nls)
    }
    
  }
  
  if (is.null(sb.ind))
  {
    print("No shape breaks")
  }
  return(sb.ind)
  
}

#plot(w,d.in, main = "slope break location")
#points(w[sb.ind], d.in[sb.ind], col = "red")

