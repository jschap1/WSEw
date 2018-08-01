# Plot predicted vs. actual z0 over many cross sections

result <- NULL
seg <- 1
while (is.null(result))
{
  result <- slope_break(xWSEw[[seg]], exposure = 1, thres = 0.015)
  seg <- seg + 1
}

# Calculate error of the estimate (absolute error w.r.t. b.min)

z0 <- vector(length = nseg)
for (seg in 1:nseg)
{
  result <- slope_break(xWSEw[[seg]], exposure = 0.6, thres = 0.015, win = 4)
  if (!is.null(result))
  {
    png(paste0(seg, ".png"))
    plot(WSE~w, data = xWSEw[[seg]])
    points(xWSEw[[seg]]$w[result$sb.ind], xWSEw[[seg]]$WSE[result$sb.ind], col = "red", pch = 19)
    dev.off()
    
    z0[seg] <- result$z0
  }
  
}
z0[z0==0] <- NA
z0_error <- z0 - b.min[!is.na(z0)]

plot(b.min[!is.na(z0)], z0[!is.na(z0)], asp = 1, xlab = "Truth", 
     ylab = "Predicted", main = "Thalweg elevation (m)", pch = 19)
abline(0,1)
