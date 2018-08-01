mc.nonlinear <- function(r, rWSEw, exposure, M = 100, sd_wse = 0.1, sd_w = 10,
                         plotf = FALSE, breaks = 0)
{
  # Monte Carlo approach to SWOT measurement uncertainty
  # M = number of MC simulations
  # plotf = flag for plots
  # breaks = 1 for "nonlinear slope break" method
  
  WSEw_full <- rWSEw[[r]]
  set.seed(704753262)
  wbf <- max(WSEw_full$w)
  observed.ind <- which(WSEw_full$w/wbf >= (1-exposure))
  WSEw_obs <- WSEw_full[observed.ind,]
  nn <- dim(WSEw_obs)[1] # number of observations
  
  # Make plots for visual analysis
  par(mfrow = c(1,1))
  png(paste0("p21_r_", r, "_expo_", exposure, "_WSEw_nl.png"))
  
  plot(WSE~w, data = WSEw_full, xlab = "width (m)", ylab = "WSE (m)",
       main = paste("Section", r), lty = 1, type = "l", lwd = 1, ylim = c(135,145))
  
  #plot(WSEw_full)
  #print(WSEw_full)
  
  z0 <- vector(length = M)
  sb.ind.debug <- vector(length = M)
  
  for (m in 1:M)
  {
    
    # corrupt observations
    e_WSE <- rnorm(nn, mean = 0, sd = sd_wse) # standard deviation of 0.25 m
    e_w <- rnorm(nn, mean = 0, sd = sd_w) # sd = 50 m, additive
    
    WSEw_corr <- WSEw_obs
    WSEw_corr$WSE <- WSEw_obs$WSE + e_WSE
    WSEw_corr$w <- WSEw_obs$w + e_w
    
    # Constrain observations to be positive:
    neg.ind <- which(WSEw_corr<0, arr.ind = TRUE) 
    WSEw_corr[neg.ind] <- 0
    
    if (plotf)
    {
      points(WSE~w, data = WSEw_corr, col = "cyan")
    }
    
    if (length(WSEw_corr$WSE)<2*5) 
    {
      print(paste("Fewer than 10 observed data points for reach", r, "exposure", exposure))
      next
    }
    
    # fit curve
    WSEw1 <- WSEw_corr
    
    if (breaks == 0)
    {
      fit <- fit_curve(WSEw1)
      z0[m] <- predict(fit, newdata = data.frame(w=0))
      if (plotf)
      {
        lines(WSEw1$w, predict(fit), col = "blue", lwd = 1)
        points(0, predict(fit, newdata = data.frame(w=0)), col = "blue", pch = 19, cex = 1)
      }
      
    } else if (breaks == 1)
    {
      tryCatch(
        {
          result <- fit_curve_b1(WSEw1)
          z0[m] <- result$z0
          sb.ind <- result$sb.ind
          sb.ind.debug[m] <- sb.ind
          fit.a <- result$fits[[1]]
          fit.b <- result$fits[[2]]
          if (plotf)
          {
            lines(WSEw1$w[1:sb.ind], predict(fit.a), col = "blue", lwd = 1)
            lines(WSEw1$w[sb.ind:length(WSEw1$w)], predict(fit.b), col = "blue", lwd = 1)
            points(WSEw1$w[sb.ind], WSEw1$WSE[sb.ind], pch = 19, col = "red")
            points(0, z0[m], col = "blue", pch = 19, cex = 1)
          }
        }, 
        error = function(e) 
        {
          print(paste("nlsLM failed to converge for reach", r, "exposure", exposure, "MC trial", m))
        }
      )
    }

    z0[z0==0] <- NA
  }
  
  dev.off()
  return(z0)
  
}

# ----------------------------------------------------------------------------

fit_curve <- function(WSEw1)
{
  
  tryCatch(
    {
      #fit <- nls(WSE ~ b0 + b1*w^b2,start = list(b0 = min(WSEw1$WSE), b1 = 0.01, b2 = 1), data=WSEw1,
                 #nls.control(maxiter = 200, minFactor = 1e-5))
      
      fit <- nlsLM(WSE ~ b0 + b1*w^b2, start = list(b0 = min(WSEw1$WSE), b1 = 0.01, b2 = 1), data=WSEw1)
      
    }, 
    error = function(e) {print("error")}
  )
  
  return(fit)
  
}

# ----------------------------------------------------------------------------

fit_curve_b1 <- function(WSEw, h = 4)
{
  # Fits a piecewise nonlinear curve with one breakpoint
  
  n <- dim(WSEw)[1] # number of data points
  # h <- 4 # minimum segment length
  
  LSE.best <- 1e6 # starting value
  for (i in (1+h):(n-h))
  {
    # for each possible breakpoint, find the LSE of the best fit
    
    fit.a <- nlsLM(WSE ~ a0 + a1*w^a2, start = c(a0 = min(WSEw$WSE), a1 = 1e-4, a2 = 1), 
                   data = WSEw[1:i,])
    
    a0 <- coef(fit.a)[1]
    a1 <- coef(fit.a)[2]
    a2 <- coef(fit.a)[3]
    xb <- WSEw$w[i]
    
    fit.b <- nlsLM(WSE ~ a0+a1*xb^a2 - b1*xb^b2 + b1*w^b2, 
                   start = c(b1 = 1e-4, b2 = 1), 
                   data = WSEw[(i):n,])
    
    LSE <- sum(resid(fit.a)^2) + sum(resid(fit.b)^2)
    
    if (LSE<LSE.best)
    {
      # if LSE is low, record the estimate
      LSE.best <- LSE
      fits <- list(fit.a, fit.b)
      b.ind <- i
    }
    
  }
  
  result <- list(z0 = as.numeric(coef(fit.a))[1], sb.ind = b.ind, fits = fits)
  return(result)
  
}
