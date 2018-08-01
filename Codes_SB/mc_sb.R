mc.sb <- function(r, rWSEw, exposure, M = 1000,
                  sd_wse = 0.25, sd_w = 0.25, thres = 0.015, win = 10)
{
  # Monte Carlo approach to SWOT measurement uncertainty
  # M = number of MC simulations

  # Set random seed
  set.seed(704753262)

  # Estimate bankfull width
  wbf <- max(rWSEw[[r]]$w)

  # Keep exposed data only
  unobserved.ind <- which(rWSEw[[r]]$w/wbf < (1-exposure))
  rWSEw[[r]] <- rWSEw[[r]][-unobserved.ind,]

  nn <- dim(rWSEw[[r]])[1] # number of observations

  # Make plots for visual analysis
  par(mfrow = c(1,1))
  png(paste0("p21_r_", r, "_expo_", exposure, "_WSEw.png"))
  plot(WSE~w, data = WSEw, xlab = "width (m)", ylab = "WSE (m)",
       main = "Width-WSE relationship", lty = 1, type = "l", lwd = 3)

  z0 <- vector(length = M)
  sb.ind.debug <- vector(length = M)
  for (m in 1:M)
  {

    # corrupt observations
    e_WSE <- rnorm(nn, mean = 0, sd = sd_wse) # standard deviation of 0.25 m
    e_w <- rnorm(nn, mean = 1, sd = sd_w) # cv = 0.25

    WSEw <- rWSEw[[r]]
    WSEw$WSE <- rWSEw[[r]]$WSE + e_WSE
    WSEw$w <- (e_w)*rWSEw[[r]]$w

    # find slope break
    sb.ind <- test_slope_break(WSEw$WSE, WSEw$w, thres = thres, window = win, m = FALSE)
    if (is.null(sb.ind))
    {
      sb.ind.debug[m] <- NA
      next
    } else
    {
      sb.ind.debug[m] <- sb.ind
      points(rWSEw[[4031]]$w[sb.ind], rWSEw[[4031]]$WSE[sb.ind], col = "red", pch = 19)
    }

    # fit line to points below the slope break
    WSEw1 <- WSEw[1:sb.ind,]
    lf1 <- lm(WSE~w, data = WSEw1)
    z0[m] <- lf1$coefficients[1]

    lines(WSEw1$w, predict(lf1), col = "blue", lwd = 0.5)

  }

  dev.off()
  z0[z0==0] <- NA
  return(z0)

}
