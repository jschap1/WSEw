# 7/26/2018
# Reach averaged slope break method
# This is supposed to be the semifinal version for the abstract. 

rm(list=ls())
setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Slope_Break")
load("processed_data_p21.rda")
source("sb_functions.R")
opar <- par()

# ------------------------------------------------------------------------
# Monte Carlo approach

mc.sb <- function(r, rWSEw, exposure, M = 1000,
                  sd_wse = 0.25, sd_w = 0.25, thres = 0.015, win = 10)
{
  # Monte Carlo approach to SWOT measurement uncertainty
  # M = number of MC simulations
  
  WSEw_full <- rWSEw[[r]]
  
  # Set random seed
  set.seed(704753262)
  
  # Estimate bankfull width
  wbf <- max(WSEw_full$w)
  
  # Keep exposed data only
  observed.ind <- which(WSEw_full$w/wbf > (1-exposure))
  WSEw_obs <- WSEw_full[observed.ind,]
  
  nn <- dim(WSEw_obs)[1] # number of observations
  
  # Make plots for visual analysis
  par(mfrow = c(1,1))
  png(paste0("p21_r_", r, "_expo_", exposure, "_WSEw.png"))
  plot(WSE~w, data = WSEw_full, xlab = "width (m)", ylab = "WSE (m)",
       main = "Width-WSE relationship", lty = 1, type = "l", lwd = 3)
  
  z0 <- vector(length = M)
  sb.ind.debug <- vector(length = M)
  for (m in 1:M)
  {
    
    # corrupt observations
    e_WSE <- rnorm(nn, mean = 0, sd = sd_wse) # standard deviation of 0.25 m
    e_w <- rnorm(nn, mean = 1, sd = sd_w) # cv = 0.25
    
    WSEw_corr <- WSEw_obs
    WSEw_corr$WSE <- WSEw_obs$WSE + e_WSE
    WSEw_corr$w <- (e_w)*WSEw_obs$w
    
    # find slope break
    sb.ind <- test_slope_break(WSEw_corr$WSE, WSEw_corr$w, thres = thres, 
                               window = win, m = FALSE)
    if (is.null(sb.ind))
    {
      sb.ind.debug[m] <- NA
      next
    } else
    {
      sb.ind.debug[m] <- sb.ind
      # points(WSEw_corr$w[sb.ind], WSEw_corr$WSE[sb.ind], col = "red", pch = 19)
    }
    
    # fit line to points below the slope break
    WSEw1 <- WSEw_corr[1:sb.ind,]
    lf1 <- lm(WSE~w, data = WSEw1)
    z0[m] <- lf1$coefficients[1]
    
    lines(WSEw1$w, predict(lf1), col = "blue", lwd = 0.5)
    
    # Choose the lowest observation for comparison
    # z0.low <- WSEw_corr[1]
    
  }
  
  dev.off()
  z0[z0==0] <- NA
  return(z0)
  
}

# Choose a reach 
r <- 4031

# Choose an exposure level
exposure <- 0.5

z0 <- mc.sb(r = 4031, rWSEw, exposure = 0.5, thres = 0.5)

# calculate some statistics
z0.true <- rWSEw[[r]]$WSE[1]
z0.bar <- mean(z0, na.rm = TRUE)
mse <- (1/(M-sum(is.na(z0))))*sum((z0-rep(z0.true, M))^2, na.rm = TRUE)
bias <- z0.bar - z0.true
variance <- mse-bias
stddev <- variance^0.5

# 
par(opar)
par(mfrow = c(2,1))
hist(z0, breaks = 40, main = "z0 estimates", col = "darkblue", xlim = c(133,140))
abline(v = z0.bar, col = "red", lwd = "2")
abline(v = z0.true, col = "red", lwd = "2")
legend(135, 70, legend = c(paste("Mean = ", round(z0.bar,1)),"True = ", round(z0.true,1)))
text(136.5, 40, paste("SD = ", round(stddev,1)))

# How does this compare to using the minimum observed value as the WSE?
r <- 4031
wbf <- max(rWSEw[[r]]$w)
exposure <- 0.9
unobserved.ind <- which(rWSEw[[r]]$w/wbf < (1-exposure))
rWSEw[[r]] <- rWSEw[[r]][-unobserved.ind,]
WSE.min <- min(rWSEw[[r]]$WSE) # implement this in the MC framework

# ------------------------------------------------------------------------
# Use MC approach for all exposure levels for one reach

r <- 4031
expo <- seq(0.05, 1, by = 0.05)
n_exp_levels <- length(expo)
bias <- vector(length = n_exp_levels)
variance <- vector(length = n_exp_levels)
for (m in 1:n_exp_levels)
{
  z0 <- mc.sb(r = 4031, rWSEw, exposure = expo[m], thres = 0.5)
  z0.true <- rWSEw[[r]]$WSE[1]
  z0.bar <- mean(z0, na.rm = TRUE)
  mse <- (1/(M-sum(is.na(z0))))*sum((z0-rep(z0.true, M))^2, na.rm = TRUE)
  bias[m] <- z0.bar - z0.true
  variance[m] <- mse-bias[m] # problem: MSE can be larger than bias/negative variance
  stddev <- variance^0.5
}

# How do bias and variance change with exposure level?
# How do bias and variance change with number observations (choose observations at random, or according to return period information)
# What is an appropriate variance to use for bathymetry priors?
# What is an appropriate bias correction for each exposure level?

# Plot bias and variance vs. exposure level
plot(expo, bias, col = "purple", pch = 19, xlab = "exposure (% wbf)", ylab = "(m)", 
     main = paste("z0 accuracy vs. exposure level for reach", r))
points(expo, variance^0.5, col = "blue", pch = 19)
legend("top", legend = c("Bias","Standard deviation"), pch = c(19,19), col = c("purple", "blue"))
# Something is wrong with the variance calculation.

# ------------------------------------------------------------------------
# Use MC approach for all reaches at 80% exposure

nr <- length(rWSEw)
bias <- vector(length = nr)
variance <- vector(length = nr)
M.small <- 50
for (r in 1:nr)
{
  z0 <- mc.sb(r = r, rWSEw, exposure = 0.8, thres = 0.5, M = M.small)
  z0.true <- rWSEw[[r]]$WSE[1]
  z0.bar <- mean(z0, na.rm = TRUE)
  mse <- (1/(M.small-sum(is.na(z0))))*sum((z0-rep(z0.true, M.small))^2, na.rm = TRUE)
  bias[r] <- z0.bar - z0.true
  variance[r] <- mse-bias[r] # problem: MSE can be larger than bias/negative variance
  stddev <- variance^0.5
}

# Plot bias and variance vs. exposure level
plot(1:nr, bias, col = "purple", pch = 19, xlab = "upstream to downstream", ylab = "(m)", 
     main = paste("z0 accuracy vs. distance for exposure level", 0.8))
points(1:nr, variance^0.5, col = "blue", pch = 19)
legend("top", legend = c("Bias","Standard deviation"), pch = c(19,19), col = c("purple", "blue"))
# Something is wrong with the variance calculation.

# ------------------------------------------------------------------------
# Finding breakpoints

library(strucchange)
plot(WSE~w, rWSEw[[r]])

brkpts <- breakpoints(WSE~w, data = rWSEw[[r]], h = 30)

points(rWSEw[[r]]$w[brkpts$breakpoints], 
     rWSEw[[r]]$WSE[brkpts$breakpoints],
     col = "red", pch = 19)

# Theory says there is likely to be one breakpoint, separating low flows from within bank flows
# If the floodplain were considered, there would be another breakpoint for out-of-bank flows
# Therefore, for SB method, I will use this method to find breakpoints in the next iteration.

# ------------------------------------------------------------------------
# Processing


r <- 4031

exposure <- 0.5 # half the banks are exposed
wbf <- max(rWSEw[[r]]$w)

unobserved.ind <- which(rWSEw[[r]]$w/wbf < (1-exposure))
rWSEw[[r]] <- rWSEw[[r]][-unobserved.ind,]
wbf <- max(rWSEw[[r]]$w)

par(opar)
plot(WSE~w, data = rWSEw[[r]], xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship")

sb.ind <- test_slope_break(rWSEw[[r]]$WSE, rWSEw[[r]]$w, thres = 0.015, window = 4, m = FALSE)
points(rWSEw[[r]]$w[sb.ind], rWSEw[[r]]$WSE[sb.ind], col = "red", pch = 19)

n <- length(rWSEw[[r]]$w) # number of observed points, at the given exposure level
WSEw1 <- rWSEw[[r]][1:sb.ind,] # points below the slope break
WSEw2 <- rWSEw[[r]][sb.ind:n,] # points above the slope break

# Fit the points below the slope break
lf1 <- lm(WSE~w, data = WSEw1)
summary(lf1)
plot(lf1) # check residuals

# Modify in order to get continuity at the breakpoint
# See notes, https://www.fs.fed.us/rm/pubs/rmrs_gtr189.pdf

# Fit the points above the slope break
w.sb <- WSEw2$w[1] # width at the breakpoint
WSEw2.shifted <- WSEw2 
WSEw2.shifted$w <- WSEw2$w - w.sb
offset1 <- coef(lf1)[1] + coef(lf1)[2]*w.sb
n2 <- length(WSEw2$w)
lf2 <- lm(WSE ~ 0 + w, offset = rep(offset1, n2), data = WSEw2.shifted)
summary(lf2)
plot(lf2) # check residuals

lines(WSEw1$w, predict(lf1))
lines(WSEw2$w, predict(lf2))

# Predict z0, estimate confidence interval
pred1 <- predict(lf1, newdata = data.frame(w=0), se.fit = TRUE, interval = "prediction") 
z0 <- pred1$fit[1]

# Assume there is measurement error - try Deming regression
dem.reg <- mcreg(WSEw1$w, WSEw1$WSE, error.ratio = 1, method.reg = "Deming") # 3 data points is not enough
z0.dem <- dem.reg@para[1]

# Try a nonlinear fit
fit.nls <- nls(WSE~ c1+c2*w^c3, start = list(c1=1,c2=2,c3=3), 
               data=WSEw1, control = nls.control(maxiter = 1000))
coef(fit.nls)
wse.pred <- predict(fit.nls)
lines(WSEw1$w, wse.pred, col = "blue")

# Try a polynomial fit
fit.poly <- lm(WSE~poly(w,3), data = WSEw1)
plot(fit.poly)
plot(WSE~w, data = WSEw1)
lines(WSEw1$w, predict(fit.poly))

predict(fit.poly, newdata = data.frame(w=0), se.fit = TRUE, interval = "prediction") 
# higher order (more than 2) makes extrapolation to z0 very bad

nls.resid <- wse.pred - WSEw1$WSE
mean(nls.resid)
hist(nls.resid)
plot(nls.resid)

z0.nls <- predict(fit.nls, newdata = data.frame(w=0)) # must calculate prediction interval manually. 
# See https://www.r-bloggers.com/predictnls-part-1-monte-carlo-simulation-confidence-intervals-for-nls-models/

# How many data points (n) will actually be available? Not that many. SWOT orbit. 

# Get error on the z0 predictions
z0.true <- 134.44

# Make pretty figure
plot(WSE~w, data = rWSEw.copy[[r]])
points(rWSEw[[r]]$w[sb.ind], rWSEw[[r]]$WSE[sb.ind], col = "red", pch = 19)
lines(WSEw1$w, predict(lf1))
lines(WSEw2$w, predict(lf2))

coefs <- coefficients(fit.nls)
z0 <- coefs[1]
err <- z0 - z0.true







