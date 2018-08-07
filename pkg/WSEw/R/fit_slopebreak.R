#' Slope Break Fit for WSE-w Relationship
#' 
#' Can do the following types of fits
#' 1. Linear with one breakpoint, continuous
#' 2. Linear with one breakpoint, noncontinuous
#' 3. Linear with several breakpoints, continous 
#' 4. Linear with several breakpoints, noncontinuous (returns just the lowest fit)
#' 5. Mersel method: linear with one breakpoint, noncontinous,  screens out "non-optimal" cross sections (returns just the lowest fit)
#' @param WSEw WSEw data (at a given level of exposure)
#' @param mersel flag for Mersel method
#' @param thres threshold for Mersel method, must be provided if using Mersel method. 0.015 has been used in the past.
#' @param thres moving average window for Mersel method, must be provided if using Mersel method. Default value is 4.
#' @param multiple_breaks one slope break or multiple slope breaks
#' @param continuity require continuity between piecewise fits?
#' @importFrom segmented seg.control segmented
#' @importFrom strucchange breakpoints breakfactor
#' @export

fit_slopebreak <- function(WSEw,  mersel = FALSE, thres = 0.015, 
                           window = 4, multiple_breaks = FALSE, continuity = TRUE)
{
  
  nn <- length(WSEw$w) # number of data points
  
  if (nn<=1)
  {
    # print("No observations")
    return(NULL)
  }
  
  if (mersel)
  {
    # Use original Mersel method
    sb.ind <- test_slope_break(WSEw$WSE, WSEw$w, thres = thres, window = window, m = FALSE)
    if (!is.null(sb.ind))
    {
      lf1 <- lm(WSE~w, data = WSEw[1:sb.ind,])
      fits <- list(lf1)
    } else
    {
      return(NULL)
    }
    
  } else
  {
    
    if (multiple_breaks)
    {
      # Use multiple slope breaks
      b <- breakpoints(WSE~w, data = WSEw, h=5)$breakpoints # h is the minimum number of points required for a section
      if (is.null(b)) 
      {
        b<-NA
        print("No breakpoints could be identified")
      }
      sb.ind <- b[1]

      if (!continuity)
      {
        # strucchange package does not force continuity at breakpoints
        b <- breakpoints(WSE~w, data = WSEw, h=5)
        lf <- lm(WSE ~ w*breakfactor(b), data = WSEw) 
        fits <- list(lf)
        # lf1 <- lm(WSE~w, data = WSEw[1:sb.ind,])
      
      } else
      {

        # Using segmented package to perform the fit, with breakpoints from strucchange
        sctrl <- seg.control(toll = 1e-04, it.max = 1, display = FALSE,
                             stop.if.error = TRUE, K = 10, quant = FALSE, last = TRUE, maxit.glm = 25, h = 1, 
                             n.boot=20, size.boot=NULL, gap=FALSE, jt=FALSE, nonParam=TRUE,
                             random=TRUE, powers=c(1,1), seed=NULL)
        lf<-lm(WSE~w, data = WSEw)
        lfsb<-segmented(lf, seg.Z= ~w, control=sctrl, psi=WSEw$w[b])
        fits <- list(lfsb)
        
      }
      
      # Multiple breakpoints
      # Do segmented regression with multiple (known) breakpoints
      # Require continuity (or not)
    } else
    {
      # Use one slope break
      b <- breakpoints(WSE~w, data = WSEw, breaks = 1, h=5)$breakpoints
      if (is.null(b)) 
      {
        b<-NA
        print("No breakpoints could be identified")
      }
      
      sb.ind <- b[1]
      WSEw1 <- WSEw[1:sb.ind,]
      WSEw2 <- WSEw[sb.ind:nn,]
      
      if (!continuity)
      {
        lf1 <- lm(WSE~w, data = WSEw1) # below slope break
        lf2 <- lm(WSE~w, data = WSEw2) # above slope break
      } else if (continuity)
      {
        
        lf1 <- lm(WSE~w, data = WSEw1) 
        a0 <- as.numeric(coef(lf1)[1])
        a1 <- as.numeric(coef(lf1)[2])
        wb <- WSEw$w[sb.ind]
        intercept <- a0+a1*wb
        lf2 <- lm(WSE ~ -1 + I(w-wb), data = WSEw2, offset = rep(intercept,dim(WSEw2)[1]))
        
      }
      fits <- list(lf1, lf2)
    }
  }
  
  attributes(fits) <- list(sb.ind = sb.ind) # adding the sb.ind as an output
  return(fits)
}

#' Slope Break Test
#' 
#' Tests whether there is an "optimal location" for the slope break method, as described by Mersel et al. (2013)
#'
#' @param WSE
#' @param w
#' @param thres threshold for maxdiff in slope for the WSE-w relationship to be considered linear
#' @param window
#' @param m flag for multiple slope breaks
#' @return If it is an optimal location, returns the slope break location. Otherwise returns NULL.
#' @export
test_slope_break <- function(WSE, w, thres = 0.015, window = 4, m = FALSE)
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
  
  s1.high <- s1[(nn-(window-1)):nn] 
  s.bar <- mean(s1.high, na.rm = TRUE)
  
  # compare to s1 at lower water surfaces to find a slope break
  for (j in 1:(nn-window-1))
  {
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

