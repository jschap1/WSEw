# 7/27/2018
# Slope break method, Monte Carlo implementation
# This is supposed to be the final version for the abstract. 

rm(list=ls())
setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Slope_Break/Codes_07272018")
load("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Slope_Break/Data/processed_data_p21.rda")
source("mc_sb.R")
opar <- par()

# ------------------------------------------------------------------------
# Monte Carlo approach
# For a cross section

# Choose a reach 
nseg <- length(xWSEw)
seg <- 1000

# Choose an exposure level
exposure <- 0.5

plot(WSE~w, xWSEw[[seg]])

z0 <- mc.sb(r = 4031, rWSEw = rWSEw, exposure = 0.8, M = 50, 
            plotf = TRUE, sd_wse = 0.1, sd_w = 10, nbreaks = 1)

z0 <- mc.sb(r = 4031, rWSEw = rWSEw, exposure = 0.8, M = 50, 
            plotf = TRUE, sd_wse = 0.1, sd_w = 10)

z0 <- mc.sb(r = 1000, rWSEw = xWSEw, exposure = 0.8, M = 50, 
            plotf = TRUE, sd_wse = 0.1, sd_w = 10, nbreaks = 2)

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
hist(z0, breaks = 7, main = "z0 estimates", col = "darkblue", xlim = c(133,140))
abline(v = z0.bar, col = "red", lwd = "2")
abline(v = z0.true, col = "red", lwd = "2")
legend(135, 6, legend = c(paste("Mean = ", round(z0.bar,1)),"True = ", round(z0.true,1)))
text(134, 5, paste("SD = ", round(stddev,1)))

# How does this compare to using the minimum observed value as the WSE?
r <- 4031
wbf <- max(rWSEw[[r]]$w)
exposure <- 0.9
unobserved.ind <- which(rWSEw[[r]]$w/wbf < (1-exposure))
rWSEw[[r]] <- rWSEw[[r]][-unobserved.ind,]
WSE.min <- min(rWSEw[[r]]$WSE) # implement this in the MC framework


# ------------------------------------------------------------------------
# Monte Carlo approach

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
hist(z0, breaks = 10, main = "z0 estimates", col = "darkblue", xlim = c(137.5,138.5))
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

r <- 2000
expo <- seq(0.05, 1, by = 0.05)
n_exp_levels <- length(expo)
bias <- vector(length = n_exp_levels)
variance <- vector(length = n_exp_levels)
for (m in 1:n_exp_levels)
{
  z0 <- mc.sb(r = r, rWSEw, exposure = expo[m], plotf = TRUE, M = 100, 
              sd_wse = 0.1, sd_w = 10)
  z0.true <- rWSEw[[r]]$WSE[1]
  z0.bar <- mean(z0, na.rm = TRUE)
  bias[m] <- z0.bar - z0.true
  variance[m] <- var(z0, na.rm = TRUE) # problem: MSE can be larger than bias/negative variance
  stddev <- variance^0.5
}

# How do bias and variance change with exposure level?
# How do bias and variance change with number observations (choose observations at random, or according to return period information)
# What is an appropriate variance to use for bathymetry priors?
# What is an appropriate bias correction for each exposure level?

# Plot bias and variance vs. exposure level
par(opar)
plot(expo, bias, col = "purple", pch = 19, xlab = "exposure (% wbf)", ylab = "(m)", 
     main = paste("z0 accuracy vs. exposure level for reach", r), ylim = c(-0.5,3.5))
points(expo, variance^0.5, col = "blue", pch = 19)
abline(h=0)
legend("topright", legend = c("Bias","Standard deviation"), pch = c(19,19), col = c("purple", "blue"))

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

# Plot bias and variance vs. reach
plot(1:nr, bias, col = "purple", pch = 19, xlab = "upstream to downstream", ylab = "(m)", 
     main = paste("z0 accuracy vs. distance for exposure level", 0.8))
points(1:nr, variance^0.5, col = "blue", pch = 19)
legend("top", legend = c("Bias","Standard deviation"), pch = c(19,19), col = c("purple", "blue"))
# Something is wrong with the variance calculation.

# ------------------------------------------------------------------------
# Use MC approach for all reaches at different exposure levels

nr <- length(rWSEw)
bias <- vector(length = nr)
variance <- vector(length = nr)
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


expo <- seq(0.05, 1, by = 0.05)
n_exp_levels <- length(expo)
bias <- vector(length = n_exp_levels)
variance <- vector(length = n_exp_levels)
for (m in 1:n_exp_levels)
{
  z0 <- mc.sb(r = r, rWSEw, exposure = expo[m], plotf = TRUE, M = 100, 
              sd_wse = 0.1, sd_w = 10)
  z0.true <- rWSEw[[r]]$WSE[1]
  z0.bar <- mean(z0, na.rm = TRUE)
  bias[m] <- z0.bar - z0.true
  variance[m] <- var(z0, na.rm = TRUE) # problem: MSE can be larger than bias/negative variance
  stddev <- variance^0.5
}

start.time <- Sys.time() # record start time
z0 <- mc.sb(r = r, rWSEw, exposure = expo[m], plotf = TRUE, M = 100, 
            sd_wse = 0.1, sd_w = 10, nbreaks = 2)
end.time <- Sys.time() # record end time
time.taken <- end.time - start.time # compute elapsed time
print(time.taken)
# It takes 46.81 seconds to run deming regression and unknown breakpoints, M=100
# It takes 38.48 seconds to run deming regression and known breakpoints, M=100
# It takes 20.14 seconds to run OLS and known breakpoints, M=100

# Estimated time to perform calculations
4031 # reach averaged cross sections
6035 # cross sections
19 # exposure levels
19*4031*20.14/3600 # 428 hours to run all reach average cross sections
# It would take about 3 weeks without changing something. That is much too long.

# I have how many hours between this afternoon and next time I work?
# Friday at 4 p.m. through Sunday at 8 a.m.
# About 40 hours
# How many cross sections can I run?
40*3600 # seconds
# 20 = seconds per cross section
40*3600/20
# 7200 cross sections
# Exposure levels
7200/19
# About 375 cross sections at each exposure level
# This could probably be sped up if I used fewer data points to characterize each cross section

# Plot bias and variance vs. exposure level
plot(expo, bias, col = "purple", pch = 19, xlab = "exposure (% wbf)", ylab = "(m)", 
     main = paste("z0 accuracy vs. exposure level for reach", r), ylim = c(-2,4))
points(expo, variance^0.5, col = "blue", pch = 19)
abline(h=0)
legend("topright", legend = c("Bias","Standard deviation"), pch = c(19,19), col = c("purple", "blue"))

# ------------------------------------------------------------------------
# Run for a subset of cross sections

run.ind <- round(seq(1,length(rWSEw), length.out = 10)) # It takes about 5 minutes per reach
nr <- length(run.ind)
bias <- array(dim = c(nr, n_exp_levels))
variance <- array(dim = c(nr, n_exp_levels))
z0.bar <- array(dim = c(nr, n_exp_levels))
z0.true <- vector(length = nr)
expo <- seq(0.05, 0.95, by = 0.05)
n_exp_levels <- length(expo)
z0 <- list(length = nr*n_exp_levels)
ind <- 1
for (r in 1:nr)
{
  for (m in 1:n_exp_levels)
  {
    z0[[ind]] <- mc.sb(r = run.ind[r], rWSEw, exposure = expo[m], plotf = TRUE, M = 100, 
                sd_wse = 0.1, sd_w = 10)
    z0.true[r] <- rWSEw[[run.ind[r]]]$WSE[1]
    z0.bar[r,m] <- mean(z0[[ind]], na.rm = TRUE)
    bias[r,m] <- z0.bar[r,m] - z0.true[r]
    variance[r,m] <- var(z0[[ind]], na.rm = TRUE) 
    ind <- ind + 1
  }
}
save(z0, z0.true, z0.bar, bias, variance, file = "reach_avg_results.rda")

########################
# Given a known cross section shape, how well can we estimate z0 at low exposure %?


########################

# ------------------------------------------------------------------------
# Run for a subset of cross sections

run.ind <- round(seq(1,length(xWSEw), length.out = 300)) # It takes about 5 minutes per reach
nr <- length(run.ind)
expo <- seq(0.05, 0.95, by = 0.05)
n_exp_levels <- length(expo)
bias <- array(dim = c(nr, n_exp_levels))
variance <- array(dim = c(nr, n_exp_levels))
z0.bar <- array(dim = c(nr, n_exp_levels))
z0.true <- vector(length = nr)
z0 <- list(length = nr*n_exp_levels)
ind <- 1
for (r in 1:nr)
{
  for (m in 1:n_exp_levels)
  {
    z0[[ind]] <- mc.sb(r = run.ind[r], rWSEw = xWSEw, exposure = expo[m], plotf = TRUE, M = 100, 
                       sd_wse = 0.1, sd_w = 10)
    z0.true[r] <- xWSEw[[run.ind[r]]]$WSE[1]
    z0.bar[r,m] <- mean(z0[[ind]], na.rm = TRUE)
    bias[r,m] <- z0.bar[r,m] - z0.true[r]
    variance[r,m] <- var(z0[[ind]], na.rm = TRUE) 
    ind <- ind + 1
  }
}
save(z0, z0.true, z0.bar, bias, variance, file = "cross_sections_results.rda")

# ------------------------------------------------------------------------
# Get a feel for how many breakpoints there are by running breakpoints on all cross sections

nb <- vector(length = length(run.ind))
for (r in 1:length(run.ind))
{
  b <- breakpoints(WSE~w, data = rWSEw[[run.ind[r]]], h=0.2)$breakpoints
  nb[r] <- length(b)
}
nb
boxplot(nb, main = "number of breakpoints")

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
# Try running mc.sb in parallel
library(parallel)
numcores <- detectCores() - 1 # don't use all the cores!
cl <- makeCluster(numcores) # initiate cluster
start.time <- Sys.time() # record start time

# zz <- lapply(X = seq(1:2), FUN = mc.sb, rWSEw = rWSEw, exposure = 0.5, plotf = TRUE, M = 10)

listinput <- list(r = 1:20, rWSEw = rWSEw, exposure = 0.5, plotf = TRUE, M = 10)

# str(rep(list(rWSEw, exposure = 0.5, plotf = TRUE, M = 50),4))

result <- parLapply(cl = cl, listinput, fun = mc.sb)

end.time <- Sys.time() # record end time
time.taken <- end.time - start.time # compute elapsed time
print(time.taken)

# Return control of resources to OS
stopCluster(cl)

# ------------------------------------------------------------------------
# Plot histograms of error
M <- 1e3
e_wse <- rnorm(n = M, mean = 0, sd = 0.1)
e_w <- rnorm(n = M, mean = 0, sd = 10)
par(mfrow = c(2,1))
hist(e_wse, breaks = 40, main = "WSE error (m)", col = "darkblue")
hist(e_w, breaks = 40, main = "w error (m)", col = "darkblue")

