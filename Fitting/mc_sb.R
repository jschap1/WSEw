mc.sb <- function(r, rWSEw, exposure, M = 100, sd_wse = 0.1, sd_w = 10,
                  plotf = FALSE, nbreaks = NULL)
{
  # Monte Carlo approach to SWOT measurement uncertainty
  # Uses strucchange::breakpoints to find slope breaks
  # M = number of MC simulations
  # plotf = flag for plots
  
  WSEw_full <- rWSEw[[r]]
  
  # Set random seed
  set.seed(704753262)
  
  # Estimate bankfull width
  wbf <- max(WSEw_full$w)
  
  # Keep exposed data only
  observed.ind <- which(WSEw_full$w/wbf > (1-exposure))
  # if (length(observed.ind)<2){next}
  WSEw_obs <- WSEw_full[observed.ind,]
  
  nn <- dim(WSEw_obs)[1] # number of observations
  
  if (plotf)
  {
    # Make plots for visual analysis
    par(mfrow = c(1,1))
    png(paste0("p21_r_", r, "_expo_", exposure, "_WSEw.png"))
    plot(WSE~w, data = WSEw_full, xlab = "width (m)", ylab = "WSE (m)",
         main = paste("Section", r), lty = 1, type = "l", lwd = 1)
  }

  z0 <- vector(length = M)
  z0.low <- vector(length = M)
  sb.ind.debug <- vector(length = M)
  
  for (m in 1:M)
  {
    
    # corrupt observations
    e_WSE <- rnorm(nn, mean = 0, sd = sd_wse) # standard deviation of 0.25 m
    # e_w <- rnorm(nn, mean = 1, sd = sd_w) # cv = 0.25, multiplicative
    e_w <- rnorm(nn, mean = 0, sd = sd_w) # sd = 50 m, additive
    
    WSEw_corr <- WSEw_obs
    WSEw_corr$WSE <- WSEw_obs$WSE + e_WSE
    # WSEw_corr$w <- (e_w)*WSEw_obs$w
    WSEw_corr$w <- WSEw_obs$w + e_w
    
    if (plotf)
    {
      points(WSE~w, data = WSEw_corr, col = "cyan")
      #lines(WSE~w, data = WSEw_corr, col = "cyan")
    }
    
    if (length(WSEw_corr$WSE)<2*5) {next}
    
    # find slope break
    if (is.null(nbreaks))
    { # let strucchange find the number of breakpoints
      b <- breakpoints(WSE~w, data = WSEw_corr, h=5)$breakpoints
    } else
    { # specify the number of breakpoints
      b <- breakpoints(WSE~w, data = WSEw_corr, breaks = nbreaks, h=5)$breakpoints
    }
    # print(b)
    sb.ind <- b[1] # take the lowest breakpoint
    # print(sb.ind)
    if (is.null(sb.ind) | is.na(sb.ind))
    {
      sb.ind.debug[m] <- NA
      next
    } else
    {
      sb.ind.debug[m] <- sb.ind
    }
    
    # fit line to points below the slope break
    WSEw1 <- WSEw_corr[1:sb.ind,]
    lf1 <- lm(WSE~w, data = WSEw1)
    z0[m] <- lf1$coefficients[1]
    
    # What if I used Deming regression?
    #dem.reg <- mcreg(WSEw1$w, WSEw1$WSE, error.ratio = sd_wse/sd_w, method.reg = "Deming") 
    #z0[m] <- dem.reg@para[1,1]
    #s <- dem.reg@para[2,1]
    
    if (plotf)
    {
      
      #lines(WSEw1$w, WSEw1$w*s+z0[m], col = "blue", lwd = 1)
      #points(0, z0[m], col = "blue", pch = 19, cex = 1)
      points(WSEw_corr$w[b], WSEw_corr$WSE[b], col = "red", pch = 19)
      lines(WSEw1$w, predict(lf1), col = "blue", lwd = 1)
      points(0, predict(lf1, newdata = data.frame(w=0)), col = "blue", pch = 19, cex = 1)
    }
    
    # Choose the lowest observation for comparison
    z0.low[m] <- WSEw_corr$WSE[1]
    # print(z0.low[m])
    
    # Display progress
    if (m%%10 == 0)
    {
      print(m)
    }
    
  }
  
  if (plotf)
  {
    dev.off()
  }
  
  print(paste("the average minimum observed WSE was", mean(z0.low, na.rm = TRUE)))
  z0[z0==0] <- NA
  
  # result <- list(z0 = z0, sb.ind = sb.ind.debug, fits = lf1, z0.low)
  
  #return(result)
  return(z0)
  
}

# ---------------------------------------------------------------------------------------------------------------- 
# How many breakpoints are there and where do they occur?

find_lin_breakpoint <- function(WSEw, nbreaks = NULL)
{
  # Can pre-specify the number of slopebreaks or not
  # Requires strucchange package::breakpoints
  
  # find slope break
  if (is.null(nbreaks))
  { # let strucchange find the number of breakpoints
    b <- breakpoints(WSE~w, data = WSEw, h=5)$breakpoints
  } else
  { # specify the number of breakpoints
    b <- breakpoints(WSE~w, data = WSEw, breaks = nbreaks, h=5)$breakpoints
  }

  sb.ind <- b[1] # take the lowest breakpoint
  
  if (is.null(b))
  {
    b <- NA
  }
  
  return(b)
  # return(sb.ind)
}



