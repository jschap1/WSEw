# Main code for fitting cross sections
# 7/27/2018
#
# Fits reach average and individual cross sections
# Uses linear, slope break, multiple slope break, nonlinear, and nonlinear slope break methods

rm(list=ls())
setwd("/Users/jschap/Desktop/Cross_Sections/")
opar <- par()
library(minpack.lm)

source("/Users/jschap/Desktop/Cross_Sections/Codes/Fitting/mc_sb.R")
source("/Users/jschap/Desktop/Cross_Sections/Codes/Fitting/mc_linear.R")
source("/Users/jschap/Desktop/Cross_Sections/Codes/Fitting/mc_nonlinear.R")

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
sbsavename <- "r_nr76_expo20_500m_0.05_sb.rda" # slope break
lsavename <- "r_nr76_expo20_500m_0.05_l.rda" # linear
nlsavename <- "r_nr76_expo20_500m_0.05_nl.rda" # nonlinear
nlsbsavename <- "xs_nr76_expo10_500m_0.05_nlsb.rda" # shape break

lmsavename <- "xs_nr3_expo10_500m_0.05_lm.rda" # linear (Mersel)
sbmsavename <- "xs_nr3_expo10_500m_0.05_sbm.rda" # slope break (Mersel)

# END INPUTS
# ----------------------------------------------------------------------------------------------------------------------------
# Do slope break method

# If there is one slope break, where does it occur?
sb.ind <- vector(length = nr)
sb.expo <- vector(length = nr)
n.obs <- dim(rWSEw[[1]])[1]
w.bf <- unlist(lapply(rWSEw, max))
for (r in 1:nr)
{
  #sb.ind[r] <- find_lin_breakpoint(rWSEw[[r]], nbreaks = 1)
  sb.ind[r] <- find_lin_breakpoint(rWSEw[[r]])
  wsb <- rWSEw[[r]][sb.ind[r],2] # width of slope break
  sb.expo[r] <- 1-wsb/w.bf[r] # bankfull width
  # at what fraction of bankfull width do the slope breaks occur?
}
hist(sb.expo, main = "SB occurrence", col = "darkblue", xlab = "channel exposure (%)")

# How does this change for individual cross sections (not reach averaged)?
n.xs <- length(xWSEw)
sb.ind <- vector(length = n.xs)
sb.expo <- vector(length = n.xs)
n.obs <- dim(xWSEw[[1]])[1]
w.bf <- unlist(lapply(xWSEw, max))
for (r in 1:n.xs)
{
  sb.ind[r] <- find_lin_breakpoint(xWSEw[[r]])
  wsb <- xWSEw[[r]][sb.ind[r],2] # width of slope break
  sb.expo[r] <- 1-wsb/w.bf[r] # bankfull width
}
hist(sb.expo, main = "SB occurrence", col = "darkblue", xlab = "channel exposure (%)")


# Perform slope break method for all cross sections
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
    z0[[ind]] <- mc.sb(r = run.ind[r], rWSEw = rWSEw, exposure = expo[m], plotf = TRUE, M = 100,# rWSEw argument
                       sd_wse = 0.1, sd_w = 10)
    z0.true[r] <- rWSEw[[run.ind[r]]]$WSE[1] # rWSEw argument
    z0.bar[r,m] <- mean(z0[[ind]], na.rm = TRUE)
    bias[r,m] <- z0.bar[r,m] - z0.true[r]
    variance[r,m] <- var(z0[[ind]], na.rm = TRUE) 
    ind <- ind + 1
    print(m)
  }
  print(paste0("Progress: ", 100*round(r/nr, 2), "%"))
}
save(z0, z0.true, z0.bar, bias, variance, file = file.path(saveloc, sbsavename))

# ----------------------------------------------------------------------------------------------------------------------------
# Do linear method

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
    z0[[ind]] <- mc.linear(r = run.ind[r], rWSEw = rWSEw, exposure = expo[m], plotf = FALSE, M = 100, # rWSEw argument
                       sd_wse = 0.1, sd_w = 10)
    z0.true[r] <- rWSEw[[run.ind[r]]]$WSE[1] # rWSEw argument
    z0.bar[r,m] <- mean(z0[[ind]], na.rm = TRUE)
    bias[r,m] <- z0.bar[r,m] - z0.true[r]
    variance[r,m] <- var(z0[[ind]], na.rm = TRUE) 
    ind <- ind + 1
  }
  print(paste0("Progress: ", 100*round(r/nr, 2), "%"))
}
save(z0, z0.true, z0.bar, bias, variance, file = file.path(saveloc, lsavename))

# ----------------------------------------------------------------------------------------------------------------------------
# Do nonlinear method

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
    z0[[ind]] <- mc.nonlinear(r = run.ind[r], rWSEw = rWSEw, exposure = expo[m], plotf = TRUE, M = 100, # rWSEw argument
                           sd_wse = 0.1, sd_w = 10, breaks = 0)
    z0.true[r] <- rWSEw[[run.ind[r]]]$WSE[1] # rWSEw argument
    z0.bar[r,m] <- mean(z0[[ind]], na.rm = TRUE)
    bias[r,m] <- z0.bar[r,m] - z0.true[r]
    variance[r,m] <- var(z0[[ind]], na.rm = TRUE) 
    ind <- ind + 1
  }
  print(paste0("Progress: ", 100*round(r/nr, 2), "%"))
}
save(z0, z0.true, z0.bar, bias, variance, file = file.path(saveloc, nlsavename))

# ----------------------------------------------------------------------------------------------------------------------------
# Do nonlinear slope break method

n_exp_levels <- length(expo) # takes about 1.5 hours for 76 cross sections with 0.05 discretization
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
    z0[[ind]] <- mc.nonlinear(r = run.ind[r], rWSEw = rWSEw, exposure = expo[m], plotf = TRUE, M = 100, # rWSEw argument
                              sd_wse = 0.1, sd_w = 10, breaks = 1)
    z0.true[r] <- rWSEw[[run.ind[r]]]$WSE[1] # rWSEw argument
    z0.bar[r,m] <- mean(z0[[ind]], na.rm = TRUE)
    bias[r,m] <- z0.bar[r,m] - z0.true[r]
    variance[r,m] <- var(z0[[ind]], na.rm = TRUE) 
    ind <- ind + 1
  }
  print(paste0("Progress: ", 100*round(r/nr, 2), "%"))
}
save(z0, z0.true, z0.bar, bias, variance, file = file.path(saveloc, nlsbsavename))


# ----------------------------------------------------------------------------------------------------------------------------
# Plot average error for each method for a given exposure, over all the cross sections

par(mfrow = c(1,1), mar = c(5,5,2,5))
df <- data.frame(variance = colMeans(variance, na.rm = TRUE), 
                 bias = colMeans(bias, na.rm = TRUE), exp = expo)

with(df, plot(100*exp, bias, main = "Bottom Elevation Prediction Accuracy",
              xlab = "Channel exposure (%)", 
              ylab = "Bias (m)", pch = 19, lwd = 2, col = "red", ylim = c(-10,2)))
abline(0,0)

par(new = TRUE)
with(df, plot(100*exp, variance, pch = 25, col = "blue", lwd = 2, 
              axes = FALSE, xlab = NA, ylab = NA))

axis(side = 4)
mtext(side = 4, line = 3, "Variance")
legend("top", col = c("red", "blue"), pch = c(19, 25),
       legend=c("Bias","Variance"))


# ----------------------------------------------------------------------------------------------------------------------------
# Plot on one set of axes

load(file.path(saveloc, lsavename))
z0l <- z0
biasl <- bias
variancel <- variance

load(file.path(saveloc, sbsavename))
z0sb <- z0
biassb <- bias
variancesb <- variance

load(file.path(saveloc, nlsavename))
z0lnl<- z0
biasnl <- bias
variancenl <- variance

load(file.path(saveloc, nlsbsavename))
z0lnlsb<- z0
biasnlsb <- bias
variancenlsb <- variance

# Bias
plot(expo, colMeans(biasl, na.rm = TRUE), 
     main = "Bias comparison", xlab = "channel exposure (%)", 
     ylab = "bias (m)", col = "blue", pch = 19, ylim = c(-2, 1))
abline(0,0)
points(expo, colMeans(biassb, na.rm = TRUE), pch = 19, col = "darkgreen")
points(expo, colMeans(biasnl, na.rm = TRUE), pch = 19, col = "orange")
points(expo, colMeans(biasnlsb, na.rm = TRUE), pch = 19, col = "black")
legend("topleft", legend = c("Linear","Slope break","Nonlinear", "NLSB"), 
       col = c("blue","darkgreen","orange","black"), pch = c(19,19,19,19), ncol = 2)
 
# Variance
plot(expo, colMeans(variancel, na.rm = TRUE), 
     main = "Variance comparison", xlab = "channel exposure (%)", 
     ylab = "Variance (m^2)", col = "blue", pch = 19)
abline(0,0)
points(expo, colMeans(variancesb, na.rm = TRUE), pch = 19, col = "darkgreen")
points(expo, colMeans(variancenl, na.rm = TRUE), pch = 19, col = "orange")
points(expo, colMeans(variancenlsb, na.rm = TRUE), pch = 19, col = "black")
legend("top", legend = c("Linear","Slope break","Nonlinear","NLSB"), 
       col = c("blue","darkgreen","orange","black"), pch = c(19,19,19,19), ncol = 2)


