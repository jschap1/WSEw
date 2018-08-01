# Brute force breakpoint selection using Levenberg Marquadt algorithm
# Jacob Schaperow
# 7/31/2018

# Get WSE-w data
# WSEw <- data.frame(w=w, WSE=WSE)
seg <- 76
WSEw <- xWSEw[[seg]]

# Specify number of breakpoints
nb <- 1

# ------------------------------------------------------------------------------------
# Zero breakpoints

nb <- 0

# Perform the fit
fit <- nlsLM(WSE ~ a0 + a1*w^a2, start = c(a0 = min(WSEw$w), a1 = 1, a2 = 1), data = WSEw)
LSE <- sum(resid(fit)^2)
b.ind <- 1 # index of breakpoint

# Plot fit
par(mfrow = c(1,1))
plot(WSE~w, data = WSEw, main = paste("segment", seg, "LSE =", round(LSE,2)))
points(WSEw$w[b.ind], WSEw$WSE[b.ind], pch = 19, col = "red")
lines(WSEw$w, predict(fit))

# ------------------------------------------------------------------------------------
# One breakpoint

n <- dim(WSEw)[1] # number of data points
h <- 4 # minimum segment length

LSE.best <- 1e6 # starting value
for (i in (1+h):(n-h))
{
  # for each possible breakpoint, find the LSE of the best fit
  
  fit.a <- nlsLM(WSE ~ a0 + a1*w^a2, start = c(a0 = min(WSEw$w), a1 = 1e-4, a2 = 1), 
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

# Plot fit
par(mfrow = c(1,1))
plot(WSE~w, data = WSEw, main = paste("segment", seg, "LSE =", round(LSE,2)))
points(WSEw$w[b.ind], WSEw$WSE[b.ind], pch = 19, col = "red")
lines(WSEw$w[1:b.ind], predict(fits[[1]]))
lines(WSEw$w[(b.ind):n], predict(fits[[2]]))


# ------------------------------------------------------------------------------------
# Two breakpoints - start with just one breakpoint, see if it helps compared to no breakpoints




