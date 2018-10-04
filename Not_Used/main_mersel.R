# Main code for fitting cross sections - using unmodified Mersel method as a baseline
# 8/1/2018
#
# Fits reach average and individual cross sections
# Uses linear and slope break methods following Mersel, 2012

# --------------------------------------------------------------------------------------------------
# Set up environment

rm(list=ls())
opar <- par()

# Ideally, these packages would be imported by the WSEw package.
# Better yet, it wouldn't require this many libraries; it is a lot.
library(rgdal) 
library(maptools)
library(raster)
library(strucchange)
library(segmented)
library(minpack.lm)
library(WSEw)

setwd("/Users/jschap/Desktop/Cross_Sections")

# --------------------------------------------------------------------------------------------------
# Load data

# Cross section geometry and WSE-w relationships
processed_data <- "/Users/jschap/Desktop/Cross_Sections/Data/Processed_Data/processed_data_p21_sl_5m_hires.rda"
load(processed_data)

# Exposure levels
expo <- seq(0.05, 0.95, by = 0.05)

nr <- length(rWSEw)

# Make observations
WSEw_obs <- observe(rWSEw[[r]], exposure = 0.5, sd_wse = 0, sd_w = 0)

# Use linear method
lf.mersel <- fit_linear(WSEw_obs, mersel = TRUE, thres = 1e5) 
lf.mersel <- fit_linear(WSEw_obs, mersel = TRUE, thres = 0) 

# Find optimal locations
thres <- seq(0, 0.15, by = 0.015)
thres <- seq(0.5, 1, by = 0.05)
thres <- 1e4
nr <- length(rWSEw)
opt.lin <- vector(length = nr)
opt.lin.frac <- vector(length = length(thres))
ind <- 1
for (t in thres)
{
  for (r in 1:nr)
  {
    WSEw_obs <- observe(rWSEw[[r]], exposure = 1, sd_wse = 0, sd_w = 0)
    opt.lin[r] <- test_linear(WSEw_obs$WSE, WSEw_obs$w, thres = t)
  }
  opt.lin.ind <- which(opt.lin)
  opt.lin.frac[ind] <- sum(opt.lin)/nr
  print(paste("progress:", round(ind/length(thres)*100,2), "%"))
  ind <- ind + 1
}
plot(thres, opt.lin.frac, ylab = "Optimal location fraction", xlab = "maxdiff threshold", 
     main = "Linear Method") # threshold, optimal fraction

# Do fits at each threshold value

n_thres_levels <- length(thres)
bias <- array(dim = c(nr, n_thres_levels))
z0.true <- vector(length = nr)
z0 <- array(dim = c(nr, n_thres_levels))
r2 <- array(dim = c(nr, n_thres_levels))
SSE <- array(dim = c(nr, n_thres_levels))
for (r in 1:nr)
{
  ind <- 1
  WSEw_obs <- observe(rWSEw[[r]], exposure = 1, sd_wse = 0, sd_w = 0) #!
  z0.true[r] <- rWSEw[[r]]$WSE[1] #!
  for (t in thres)
  {
    lf <- fit_linear(WSEw_obs, mersel = TRUE, thres = t)
    if (!is.null(lf))
    {
      z0[r,ind] <- coef(lf[1])[1] # reach, threshold
      r2[r,ind] <- summary(lf)$r.squared
      SSE[r, ind] <- sum(lf$residuals^2)
      bias[r,ind] <- z0[r,ind] - z0.true[r]
    }
    ind <- ind + 1
  }
  print(paste("reach", r, "of", nr))
}


# ----------------------------------------------------------------------------------------------------------------------------
# Plot results for linear method

par(mfrow = c(1,1), mar = c(5,5,2,5))
df <- data.frame(thres = thres, 
                 opt = opt.lin.frac, bias = colMeans(bias, na.rm = TRUE))
with(df, plot(thres, 100*opt, main = "Linear method",
              xlab = "maxdiff threshold", 
              ylab = "Optimum locations (%)", pch = 19, lwd = 2, col = "red"))
par(new = TRUE)
with(df, plot(thres, bias, pch = 25, col = "blue", lwd = 2, 
              axes = FALSE, xlab = NA, ylab = NA))
axis(side = 4)
mtext(side = 4, line = 3, "Bias (m)")
legend("top", col = c("red", "blue"), pch = c(19, 25),
       legend=c("Optimum locations","Bias"))

# Suggests a maxdiff threshold of 0.6

r2.bar <- colMeans(r2, na.rm = TRUE)
SSE.bar <- colMeans(SSE, na.rm = TRUE)

par(mfrow = c(3,1))
plot(thres, opt.lin.frac, main = "maxdiff", type="l")
plot(thres, r2.bar, main = "R2", type="l")
plot(thres, SSE.bar, main = "SSE", type="l")

hist(r2)
hist(SSE)

quants <- quantile(r2, c(1/3, 2/3))

# Plot reach-average cross section and fit for a "good" fit
inds <- which(r2>quants[2])
WSEw <- rWSEw[[inds[1]]]
plot(WSE~w, data = WSEw, main = "Relatively good fit")
lf <- fit_linear(WSEw, mersel = TRUE, thres = 1)
lines(WSEw$w, fitted(lf))
bias <- coef(lf)[1] - z0.true[inds[1]]
legend("topleft", legend = c(paste("z0 bias =", round(bias,2))))

# "OK" fit
inds <- which(r2>quants[1] & r2<quants[2])
WSEw <- rWSEw[[inds[1]]]
plot(WSE~w, data = WSEw, main = "OK fit")
lf <- fit_linear(WSEw, mersel = TRUE, thres = 1)
lines(WSEw$w, fitted(lf))
bias <- coef(lf)[1] - z0.true[inds[1]]
legend("topleft", legend = c(paste("z0 bias =", round(bias,2))))

# "Bad" fit
inds <- which(r2<quants[1])
WSEw <- rWSEw[[inds[1]]]
plot(WSE~w, data = WSEw, main = "Relatively bad fit")
lf <- fit_linear(WSEw, mersel = TRUE, thres = 1)
lines(WSEw$w, fitted(lf))
bias <- coef(lf)[1] - z0.true[inds[1]]
legend("topleft", legend = c(paste("z0 bias =", round(bias,2))))

# ----------------------------------------------------------------------------------------------------------------------------


plot(WSE~w, data = xWSEw[[10]])
# Basically, there are always slope breaks for the reach averaged cross sections, so the linear method does not detect 
# optimal locations until the maxdiff threshold is increased to something like 0.8, which is fairly high.
# The story may be different for the individual cross section scale.

sbf <- fit_slopebreak(WSEw_obs, mersel = TRUE, thres = 0.015, window = 4)

# Reaches
# choose a number of reaches or all reaches
# run.ind <- round(seq(1,length(xWSEw), length.out = 3))
run.ind <- 1:length(rWSEw)
nr <- length(run.ind)

# Save names
saveloc <- "/Users/jschap/Desktop/Cross_Sections/Data/Fitting_Results/"
lmsavename <- "xs_nr76_expo20_500m_0.05_lm.rda" # linear (Mersel)
sbmsavename <- "xs_nr76_expo20_500m_0.05_sbm.rda" # slope break (Mersel)

# END INPUTS
# ----------------------------------------------------------------------------------------------------------------------------
# Do slope break method

# ----------------------------------------------------------------------------------------------------------------------------
# Do linear method

n_exp_levels <- length(expo)
bias <- array(dim = c(nr, n_exp_levels))
variance <- array(dim = c(nr, n_exp_levels))
z0.bar <- array(dim = c(nr, n_exp_levels))
z0.true <- vector(length = nr)
z0 <- list(length = nr*n_exp_levels)

optimal_locations <- vector(length = n_exp_levels) # optimal locations with no error

m <- 19 # full exposure
ind <- 1

for (r in 1:nr)
{
  z0[[ind]] <- linear_mersel(r = run.ind[r], rWSEw = rWSEw, 
                             exposure = expo[m], M = 1, sd_wse = 0, sd_w = 0, thres = 0.03)
  z0.true[r] <- rWSEw[[run.ind[r]]]$WSE[1]
  z0.bar[r,m] <- mean(z0[[ind]], na.rm = TRUE)
  bias[r,m] <- z0.bar[r,m] - z0.true[r]
  variance[r,m] <- var(z0[[ind]], na.rm = TRUE) 
  ind <- ind + 1
}

which(!is.na(unlist(z0))) # index of optimal reaches
optimal_locations[m] <- sum(!is.na(unlist(z0)))
print(optimal_locations)

# Plot the optimal locations
r.opt <- 57
plot(WSE~w, rWSEw[[r.opt]], main = paste("Reach", r.opt))

# save(z0, z0.true, z0.bar, bias, variance, file = file.path(saveloc, lmsavename))

# How does the MAE optimality criterion compare with something different, like SSE or R2?
# What are the R2 values for these "optimal" reaches?

lf1 <- lm(WSE~w, rWSEw[[r.opt]])
summary(lf1)
lines(rWSEw[[r.opt]]$w, predict(lf1), col = "blue")
SSE <- sum(lf1$residuals^2)

# Calculate SSE, R2, and Mersel criterion for each reach and tabulate them
SSE <- vector(length = nr)
R2 <- vector(length = nr)
max.diff <- vector(length = nr)
for (r in 1:nr)
{
  lf1 <- lm(WSE~w, rWSEw[[r]])
  SSE[r] <- sum(lf1$residuals^2)
  R2[r] <- summary(lf1)$r.squared
  max.diff[r] <- get_maxdiff(rWSEw[[r]]$WSE, rWSEw[[r]]$w)
}
df <- data.frame(SSE = SSE, R2 = R2, maxdiff = max.diff)
plot(df)
# Note: maxdiff is likely a good criterion because it prevents choosing a WSE-w relationship
# that is mostly linear, but becomes nonlinear at some point
#
# However, it would be nice to have some intuition on which value of maxdiff to specify. 
# For example, choosing maxdiff = 0.05 tends to keep WSE-w fits with SSE<4 and R2>0.94
# Mersel chose a value of 0.015 arbitrarily.

#### --------------------------------------------------------------------------------
# For updating the WSEw package

library(devtools)
library(roxygen2)
setwd("Codes/pkg")
rm(list=ls())
getwd()
setwd("WSEw")
setwd("..")
document()
install("WSEw")
library(WSEw)
?resample_polyline
?auto_transects
resample_polyline
