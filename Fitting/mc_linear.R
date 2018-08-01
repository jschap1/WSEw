mc.linear <- function(r, rWSEw, exposure, M = 100, sd_wse = 0.1, sd_w = 10,
                  plotf = FALSE)
{
  # Monte Carlo approach to SWOT measurement uncertainty
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
    png(paste0("p21_r_", r, "_expo_", exposure, "_WSEw_l.png"))
    plot(WSE~w, data = WSEw_full, xlab = "width (m)", ylab = "WSE (m)",
         main = paste("Section", r), lty = 1, type = "l", lwd = 1)
  }
  
  z0 <- vector(length = M)

  for (m in 1:M)
  {
    
    # corrupt observations
    e_WSE <- rnorm(nn, mean = 0, sd = sd_wse) # standard deviation of 0.25 m
    e_w <- rnorm(nn, mean = 0, sd = sd_w) # sd = 50 m, additive
    
    WSEw_corr <- WSEw_obs
    WSEw_corr$WSE <- WSEw_obs$WSE + e_WSE
    WSEw_corr$w <- WSEw_obs$w + e_w
    
    if (plotf)
    {
      points(WSE~w, data = WSEw_corr, col = "cyan")
    }
    
    if (length(WSEw_corr$WSE)<2*5) {next}
    
    # fit line
    WSEw1 <- WSEw_corr
    lf1 <- lm(WSE~w, data = WSEw1)
    z0[m] <- lf1$coefficients[1]
    
    if (plotf)
    {
      lines(WSEw1$w, predict(lf1), col = "blue", lwd = 1)
      points(0, predict(lf1, newdata = data.frame(w=0)), col = "blue", pch = 19, cex = 1)
    }
    
  }
  
  if (plotf)
  {
    dev.off()
  }
  
  z0[z0==0] <- NA
  return(z0)
  
}

# ----------------------------------------------------------------------------------------------------------------------------

linear_mersel <- function(r, rWSEw, exposure, M = 100, sd_wse = 0.1, sd_w = 10, thres = 0.015)
{
  # Uses the exact "linear" method of Mersel et al. (2013)
  # Only difference is that measurement error is included in an MC simulation framework

  WSEw_full <- rWSEw[[r]]
  
  # Set random seed
  set.seed(704753262)
  
  # Estimate bankfull width
  wbf <- max(WSEw_full$w)
  
  # Keep exposed data only
  observed.ind <- which(WSEw_full$w/wbf > (1-exposure))
  WSEw_obs <- WSEw_full[observed.ind,]
  
  n.obs <- dim(WSEw_obs)[1]

  z0 <- vector(length = M)
  
  for (m in 1:M)
  {
    
    # corrupt observations
    e_WSE <- rnorm(n.obs, mean = 0, sd = sd_wse) # standard deviation of 0.25 m
    e_w <- rnorm(n.obs, mean = 0, sd = sd_w) # sd = 50 m, additive
    
    WSEw_corr <- WSEw_obs
    WSEw_corr$WSE <- WSEw_obs$WSE + e_WSE
    WSEw_corr$w <- WSEw_obs$w + e_w

    if (length(WSEw_corr$WSE)<2*5) {next}
    
    # Check if an optimal slope break location
    WSEw1 <- WSEw_corr
    is.optimal <- test_linear(w = WSEw1$w, WSE = WSEw1$WSE, thres = thres)
    
    # fit line
    if (!is.null(is.optimal))
    {
      lf1 <- lm(WSE~w, data = WSEw1)
      z0[m] <- lf1$coefficients[1]
    }
    
  }
  z0[z0==0] <- NA
  return(z0)
}


