# Main code for fitting cross sections - using unmodified Mersel method as a baseline for comparison
# 8/1/2018
#
# Fits reach average and individual cross sections
# Uses linear and slope break methods following Mersel, 2012

rm(list=ls())
setwd("/Users/jschap/Desktop/Cross_Sections/")
opar <- par()

source("/Users/jschap/Desktop/Cross_Sections/Codes/Fitting/mc_sb.R")
source("/Users/jschap/Desktop/Cross_Sections/Codes/Fitting/mc_linear.R")
source("/Users/jschap/Desktop/Cross_Sections/Codes/Fitting/sb_functions.R")

# INPUTS

# Load cross section geometry and WSE-w relationships
processed_data <- "/Users/jschap/Desktop/Cross_Sections/Data/Processed_Data/processed_data_p21_500m.rda"
load(processed_data)

# Exposure levels
expo <- seq(0.05, 0.95, by = 0.05)

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


