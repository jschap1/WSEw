
# ------------------------------------------------------------------------------------------------
# Check the predicted values and models to make sure they're OK

r <- 1
k <- 16
z0.true.ra[r]
summary(z0.l[r,k,])
summary(z0.sb[r,k,])
summary(z0.nl[r,k,])
summary(z0.nlsb[r,k,])

A0.true.ra[r,k]
summary(A0.l[r,k,])
summary(A0.sb[r,k,])
summary(A0.nl[r,k,])
summary(A0.nlsb[r,k,])

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ---------------------------------------- Plotting ----------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Consider switching to Matlab and using tools from CEE251D DA project - nah, better to stick to one language

par(mfrow = c(1,1), mar = c(5,5,5,5))
par(mfrow = c(2,2))

# ------------------------------------------------------------------------------------------------
# Test for Gaussianity

# Batch test
r <- 1
k <- 19
z <- apply(s0.nlsb, c(1,2), is_normal)
z
hist(s0.nlsb[r, k,], "fd")

# The linear model z0 predictions appear Gaussian, but not the others, use quantiles in plots.

# ------------------------------------------------------------------------------------------------
# Estimate moments of sampling distributions of predicted values

# z0
z0.l.lower <- array(dim = c(nr, n_exp_levels))
z0.l.med <- array(dim = c(nr, n_exp_levels))
z0.l.upper <- array(dim = c(nr, n_exp_levels))

z0.sb.lower <- array(dim = c(nr, n_exp_levels))
z0.sb.med <- array(dim = c(nr, n_exp_levels))
z0.sb.upper <- array(dim = c(nr, n_exp_levels))

z0.nl.lower <- array(dim = c(nr, n_exp_levels))
z0.nl.med <- array(dim = c(nr, n_exp_levels))
z0.nl.upper <- array(dim = c(nr, n_exp_levels))

z0.nlsb.lower <- array(dim = c(nr, n_exp_levels))
z0.nlsb.med <- array(dim = c(nr, n_exp_levels))
z0.nlsb.upper <- array(dim = c(nr, n_exp_levels))

for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    
    q1 <- quantile(z0.l[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    z0.l.lower[r,k] <- as.numeric(q1[1])
    z0.l.med[r,k] <- as.numeric(q1[2])
    z0.l.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(z0.sb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    z0.sb.lower[r,k] <- as.numeric(q1[1])
    z0.sb.med[r,k] <- as.numeric(q1[2])
    z0.sb.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(z0.nl[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    z0.nl.lower[r,k] <- as.numeric(q1[1])
    z0.nl.med[r,k] <- as.numeric(q1[2])
    z0.nl.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(z0.nlsb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    z0.nlsb.lower[r,k] <- as.numeric(q1[1])
    z0.nlsb.med[r,k] <- as.numeric(q1[2])
    z0.nlsb.upper[r,k] <- as.numeric(q1[3])
  }
}

# A0
A0.l.lower <- array(dim = c(nr, n_exp_levels))
A0.l.med <- array(dim = c(nr, n_exp_levels))
A0.l.upper <- array(dim = c(nr, n_exp_levels))

A0.sb.lower <- array(dim = c(nr, n_exp_levels))
A0.sb.med <- array(dim = c(nr, n_exp_levels))
A0.sb.upper <- array(dim = c(nr, n_exp_levels))

A0.nl.lower <- array(dim = c(nr, n_exp_levels))
A0.nl.med <- array(dim = c(nr, n_exp_levels))
A0.nl.upper <- array(dim = c(nr, n_exp_levels))

A0.nlsb.lower <- array(dim = c(nr, n_exp_levels))
A0.nlsb.med <- array(dim = c(nr, n_exp_levels))
A0.nlsb.upper <- array(dim = c(nr, n_exp_levels))

for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    q1 <- quantile(A0.l[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A0.l.lower[r,k] <- as.numeric(q1[1])
    A0.l.med[r,k] <- as.numeric(q1[2])
    A0.l.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A0.sb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A0.sb.lower[r,k] <- as.numeric(q1[1])
    A0.sb.med[r,k] <- as.numeric(q1[2])
    A0.sb.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A0.nl[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A0.nl.lower[r,k] <- as.numeric(q1[1])
    A0.nl.med[r,k] <- as.numeric(q1[2])
    A0.nl.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A0.nlsb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A0.nlsb.lower[r,k] <- as.numeric(q1[1])
    A0.nlsb.med[r,k] <- as.numeric(q1[2])
    A0.nlsb.upper[r,k] <- as.numeric(q1[3])
  }
}

# s0
s0.l.lower <- array(dim = c(nr, n_exp_levels))
s0.l.med <- array(dim = c(nr, n_exp_levels))
s0.l.upper <- array(dim = c(nr, n_exp_levels))

s0.sb.lower <- array(dim = c(nr, n_exp_levels))
s0.sb.med <- array(dim = c(nr, n_exp_levels))
s0.sb.upper <- array(dim = c(nr, n_exp_levels))

s0.nl.lower <- array(dim = c(nr, n_exp_levels))
s0.nl.med <- array(dim = c(nr, n_exp_levels))
s0.nl.upper <- array(dim = c(nr, n_exp_levels))

s0.nlsb.lower <- array(dim = c(nr, n_exp_levels))
s0.nlsb.med <- array(dim = c(nr, n_exp_levels))
s0.nlsb.upper <- array(dim = c(nr, n_exp_levels))

for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    q1 <- quantile(s0.l[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    s0.l.lower[r,k] <- as.numeric(q1[1])
    s0.l.med[r,k] <- as.numeric(q1[2])
    s0.l.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(s0.sb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    s0.sb.lower[r,k] <- as.numeric(q1[1])
    s0.sb.med[r,k] <- as.numeric(q1[2])
    s0.sb.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(s0.nl[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    s0.nl.lower[r,k] <- as.numeric(q1[1])
    s0.nl.med[r,k] <- as.numeric(q1[2])
    s0.nl.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(s0.nlsb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    s0.nlsb.lower[r,k] <- as.numeric(q1[1])
    s0.nlsb.med[r,k] <- as.numeric(q1[2])
    s0.nlsb.upper[r,k] <- as.numeric(q1[3])
  }
}

# WP
WP.l.lower <- array(dim = c(nr, n_exp_levels))
WP.l.med <- array(dim = c(nr, n_exp_levels))
WP.l.upper <- array(dim = c(nr, n_exp_levels))

WP.sb.lower <- array(dim = c(nr, n_exp_levels))
WP.sb.med <- array(dim = c(nr, n_exp_levels))
WP.sb.upper <- array(dim = c(nr, n_exp_levels))

WP.nl.lower <- array(dim = c(nr, n_exp_levels))
WP.nl.med <- array(dim = c(nr, n_exp_levels))
WP.nl.upper <- array(dim = c(nr, n_exp_levels))

WP.nlsb.lower <- array(dim = c(nr, n_exp_levels))
WP.nlsb.med <- array(dim = c(nr, n_exp_levels))
WP.nlsb.upper <- array(dim = c(nr, n_exp_levels))

for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    q1 <- quantile(WP.l[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    WP.l.lower[r,k] <- as.numeric(q1[1])
    WP.l.med[r,k] <- as.numeric(q1[2])
    WP.l.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(WP.sb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    WP.sb.lower[r,k] <- as.numeric(q1[1])
    WP.sb.med[r,k] <- as.numeric(q1[2])
    WP.sb.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(WP.nl[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    WP.nl.lower[r,k] <- as.numeric(q1[1])
    WP.nl.med[r,k] <- as.numeric(q1[2])
    WP.nl.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(WP.nlsb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    WP.nlsb.lower[r,k] <- as.numeric(q1[1])
    WP.nlsb.med[r,k] <- as.numeric(q1[2])
    WP.nlsb.upper[r,k] <- as.numeric(q1[3])
  }
}

# A
A.l.lower <- array(dim = c(nr, n_exp_levels))
A.l.med <- array(dim = c(nr, n_exp_levels))
A.l.upper <- array(dim = c(nr, n_exp_levels))

A.sb.lower <- array(dim = c(nr, n_exp_levels))
A.sb.med <- array(dim = c(nr, n_exp_levels))
A.sb.upper <- array(dim = c(nr, n_exp_levels))

A.nl.lower <- array(dim = c(nr, n_exp_levels))
A.nl.med <- array(dim = c(nr, n_exp_levels))
A.nl.upper <- array(dim = c(nr, n_exp_levels))

A.nlsb.lower <- array(dim = c(nr, n_exp_levels))
A.nlsb.med <- array(dim = c(nr, n_exp_levels))
A.nlsb.upper <- array(dim = c(nr, n_exp_levels))

for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    q1 <- quantile(A.l[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A.l.lower[r,k] <- as.numeric(q1[1])
    A.l.med[r,k] <- as.numeric(q1[2])
    A.l.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A.sb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A.sb.lower[r,k] <- as.numeric(q1[1])
    A.sb.med[r,k] <- as.numeric(q1[2])
    A.sb.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A.nl[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A.nl.lower[r,k] <- as.numeric(q1[1])
    A.nl.med[r,k] <- as.numeric(q1[2])
    A.nl.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A.nlsb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A.nlsb.lower[r,k] <- as.numeric(q1[1])
    A.nlsb.med[r,k] <- as.numeric(q1[2])
    A.nlsb.upper[r,k] <- as.numeric(q1[3])
  }
}

# ------------------------------------------------------------------------------------------------
# Make plots of median z0, A0 along the river

par(mfrow=c(2,2))
k <- 16 # exposure level
l.pred <- data.frame(r = 1:nr, 
                     z0.med = z0.l.med[,k], z0.lower = z0.l.lower[,k], z0.upper = z0.l.upper[,k],
                     A0.med = A0.l.med[,k], A0.lower = A0.l.lower[,k], A0.upper = A0.l.upper[,k],
                     s0.med = s0.l.med[,k], s0.lower = s0.l.lower[,k], s0.upper = s0.l.upper[,k],
                     WP.med = WP.l.med[,k], WP.lower = WP.l.lower[,k], WP.upper = WP.l.upper[,k],
                     A.med = A.l.med[,k], A.lower = A.l.lower[,k], A.upper = A.l.upper[,k])
sb.pred <- data.frame(r = 1:nr, 
                      z0.med = z0.sb.med[,k], z0.lower = z0.sb.lower[,k], z0.upper = z0.sb.upper[,k],
                      A0.med = A0.sb.med[,k], A0.lower = A0.sb.lower[,k], A0.upper = A0.sb.upper[,k],
                      s0.med = s0.sb.med[,k], s0.lower = s0.sb.lower[,k], s0.upper = s0.sb.upper[,k],
                      WP.med = WP.sb.med[,k], WP.lower = WP.sb.lower[,k], WP.upper = WP.sb.upper[,k],
                      A.med = A.sb.med[,k], A.lower = A.sb.lower[,k], A.upper = A.sb.upper[,k])
nl.pred <- data.frame(r = 1:nr, 
                      z0.med = z0.nl.med[,k], z0.lower = z0.nl.lower[,k], z0.upper = z0.nl.upper[,k],
                      A0.med = A0.nl.med[,k], A0.lower = A0.nl.lower[,k], A0.upper = A0.nl.upper[,k],
                      s0.med = s0.nl.med[,k], s0.lower = s0.nl.lower[,k], s0.upper = s0.nl.upper[,k],
                      WP.med = WP.nl.med[,k], WP.lower = WP.nl.lower[,k], WP.upper = WP.nl.upper[,k],
                      A.med = A.nl.med[,k], A.lower = A.nl.lower[,k], A.upper = A.nl.upper[,k])
nlsb.pred <- data.frame(r = 1:nr, 
                        z0.med = z0.nlsb.med[,k], z0.lower = z0.nlsb.lower[,k], z0.upper = z0.nlsb.upper[,k],
                        A0.med = A0.nlsb.med[,k], A0.lower = A0.nlsb.lower[,k], A0.upper = A0.nlsb.upper[,k],
                        s0.med = s0.nlsb.med[,k], s0.lower = s0.nlsb.lower[,k], s0.upper = s0.nlsb.upper[,k],
                        WP.med = WP.nlsb.med[,k], WP.lower = WP.nlsb.lower[,k], WP.upper = WP.nlsb.upper[,k],
                        A.med = A.nlsb.med[,k], A.lower = A.nlsb.lower[,k], A.upper = A.nlsb.upper[,k])

# Do reach-averaging
l.pred.ra <- lapply(l.pred, ra, n = 2000)
sb.pred.ra <- lapply(sb.pred, ra, n = 2000)
nl.pred.ra <- lapply(nl.pred, ra, n = 2000)
nlsb.pred.ra <- lapply(nlsb.pred, ra, n = 2000)
save(l.pred.ra, sb.pred.ra, nl.pred.ra, nlsb.pred.ra, file = file.path(exp_dir, "pred_ra_80.rda")) # save predictions (ra)

# Somewhere, this is messing up the nonlinear predictions


load(file.path(exp_dir, "pred_ra_20.rda"))
load(file.path(exp_dir, "pred_ra_40.rda"))
load(file.path(exp_dir, "pred_ra_60.rda"))
load(file.path(exp_dir, "pred_ra_80.rda"))

# z0
plot(z0.true.ra, main = paste("Median z0 (m) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(100,150), lwd = 1.5,
     xlab = "reach-average cross section", ylab = "z0")
lines(l.pred.ra$z0.med, col = "red")
lines(sb.pred.ra$z0.med, col = "orange")
lines(nl.pred.ra$z0.med, col = "green")
lines(nlsb.pred.ra$z0.med, col = "blue")
legend("bottomright", legend = c("True", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A0
plot(A0.true.ra[,k], main = paste("Median A0 at", expo[k]*100,"% exposure (sq. m)"), 
     type = "l", ylim = c(0,4500), lwd = 1.5,
     xlab = "Cross section", ylab = "A0")
lines(l.pred.ra$A0.med, col = "red")
lines(sb.pred.ra$A0.med, col = "orange")
lines(nl.pred.ra$A0.med, col = "green")
lines(nlsb.pred.ra$A0.med, col = "blue")
legend("bottomright", legend = c("True", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# s0
plot(s0.true.ra*1e5, main = paste("s0 (cm/km) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(-2000,2000), 
     xlab = "reach-average cross section", ylab = "bed slope")
lines(1e5*l.pred$s0.med, col = "red")
lines(1e5*sb.pred$s0.med, col = "orange")
lines(1e5*nl.pred$s0.med, col = "green")
lines(1e5*nlsb.pred$s0.med, col = "blue")
legend("topleft", legend = c("True", "Linear","SB","NL","NLSB"), 
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

K <- 1000 # Smooth slope using a k-point moving average
s0.true.ra.smooth <- filter(s0.true.ra, sides = 2, filter = rep(1/K,K))
s0.l.smooth <- filter(l.pred$s0.med, sides = 2, filter = rep(1/K,K))
s0.sb.smooth <- filter(sb.pred$s0.med, sides = 2, filter = rep(1/K,K))
s0.nl.smooth <- filter(nl.pred$s0.med, sides = 2, filter = rep(1/K,K))
s0.nlsb.smooth <- filter(nlsb.pred$s0.med, sides = 2, filter = rep(1/K,K))

# s0 smoothed
plot(s0.true.ra.smooth*1e5, main = paste("s0 (cm/km) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(-500,500), 
     xlab = "reach-average cross section", ylab = "smoothed bed slope")
lines(1e5*s0.l.smooth, col = "red")
lines(1e5*s0.sb.smooth, col = "orange")
lines(1e5*s0.nl.smooth, col = "green")
lines(1e5*s0.nlsb.smooth, col = "blue")
legend("bottomright", legend = c("True", "Linear","SB","NL","NLSB"), 
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A
plot(A.true.ra[1:10], main = paste("Flow area (m^2) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(0,6000))
lines(l.pred$A.med, col = "red")
lines(sb.pred$A.med, col = "orange")
lines(nl.pred$A.med, col = "green")
lines(nlsb.pred$A.med, col = "blue")
legend("bottomleft", legend = c("True", "Linear","SB","NL","NLSB"), 
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# WP
plot(WP.true.ra[1:10], main = paste("Wetted perimeter (m) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(450,850))
lines(l.pred$WP.med, col = "red")
lines(sb.pred$WP.med, col = "orange")
lines(nl.pred$WP.med, col = "green")
lines(nlsb.pred$WP.med, col = "blue")
legend("topleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# ------------------------------------------------------------------------------------------------
# Make plots of average z0, A, WP, s0, A0 error at each exposure level

# z0
z0.bias <- plot_bias(expo, z0_med, z0.true.ra[1:nr], na.rm = TRUE,
                     main = "z0 bias vs. exposure level, no meas. error", ylab = "Bias (m)", ylim = c(-25,3))
legend("bottomright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# s0
s0.bias <- plot_bias(expo, s0_med, s0.true.ra[1:(nr-1)], na.rm = TRUE,
                     main = "s0 bias vs. exposure level, no meas. error", ylab = "Bias (cm/km)")
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A
A.bias <- plot_bias(expo, A_med, A.true.ra[1:nr], na.rm = TRUE,
                    main = "A bias vs. exposure level, no meas. error", ylab = "Bias (m^2)")
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# WP
WP.bias <- plot_bias(expo, WP_med, WP.true.ra[1:nr], na.rm = TRUE,
                     main = "WP bias vs. exposure level, no meas. error", ylab = "Bias (m)")
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A0
A0.bias <- array(dim = c(n_exp_levels, 4)) # cannot use plot_bias because A0 changes with exposure level
na.rm = TRUE
for (k in 1:n_exp_levels)
{
  A0.bias[k,1] <- mean(A0.l.med - A0.true.ra[,k], na.rm = na.rm)
  A0.bias[k,2] <- mean(A0.sb.med - A0.true.ra[,k], na.rm = na.rm)
  A0.bias[k,3] <- mean(A0.nl.med - A0.true.ra[,k], na.rm = na.rm)
  A0.bias[k,4] <- mean(A0.nlsb.med - A0.true.ra[,k], na.rm = na.rm)
}
A0.bias <- as.data.frame(A0.bias)
names(A0.bias) <- c("l","sb","nl","nlsb")

# Plot A0.bias vs. exposure level
par(opar)
plot(100*expo, abs(A0.bias$l), col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "Abs(A0.bias) using median as estimate (m2)", ylab = "A0.bias (m)", ylim = c(0,2000))
lines(100*expo, abs(A0.bias$sb), col = "orange")
lines(100*expo, abs(A0.bias$nl), col = "green")
lines(100*expo, abs(A0.bias$nlsb), col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

save(z0.bias, s0.bias, A.bias, WP.bias, A0.bias, file = file.path(exp_dir, "bias.rda"))

# ------------------------------------------------------------------------------------------------
# Compute RMSE

z0.rmse <- compute_rmse(z0_med, z0.true.ra)
A.rmse <- compute_rmse(A_med, A.true.ra)
WP.rmse <- compute_rmse(WP_med, WP.true.ra)
s0.rmse <- compute_rmse(s0_med, s0.true.ra)

A0.rmse <- array(dim = c(n_exp_levels, 4)) # cannot use compute_rmse because A0 changes with exposure level
na.rm = TRUE
for (k in 1:n_exp_levels)
{
  A0.rmse[k,1] <- ((1/length(A0.true.ra[,k]))*t(A0.l.med[,k] - A0.true.ra[,k])%*%(A0.l.med[,k] - A0.true.ra[,k]))^0.5
  A0.rmse[k,2] <- ((1/length(A0.true.ra[,k]))*t(A0.sb.med[,k] - A0.true.ra[,k])%*%(A0.sb.med[,k] - A0.true.ra[,k]))^0.5
  A0.rmse[k,3] <- ((1/length(A0.true.ra[,k]))*t(A0.nl.med[,k] - A0.true.ra[,k])%*%(A0.nl.med[,k] - A0.true.ra[,k]))^0.5
  A0.rmse[k,4] <- ((1/length(A0.true.ra[,k]))*t(A0.nlsb.med[,k] - A0.true.ra[,k])%*%(A0.nlsb.med[,k] - A0.true.ra[,k]))^0.5
}
A0.rmse <- as.data.frame(A0.rmse)
names(A0.rmse) <- c("l","sb","nl","nlsb")

save(z0.rmse, s0.rmse, A.rmse, WP.rmse, A0.rmse, file = file.path(exp_dir, "rmse.rda"))

# ------------------------------------------------------------------------------------------------
# Make plots of z0, A, WP, s0, A0 RMSE at each exposure level

# z0
plot(100*expo, z0.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "z0 RMSE (m)", ylab = "z0.rmse (m)",
     ylim = c(0,10))
lines(100*expo, z0.rmse$sb, col = "orange")
lines(100*expo, z0.rmse$nl, col = "green")
lines(100*expo, z0.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# s0
plot(100*expo, s0.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "s0 RMSE (cm/km)", ylab = "s0.rmse (m)")
lines(100*expo, s0.rmse$sb, col = "orange")
lines(100*expo, s0.rmse$nl, col = "green")
lines(100*expo, s0.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A
plot(100*expo, A.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "A RMSE (m2)", ylab = "A.rmse (m)",
     ylim = c(0,7000))
lines(100*expo, A.rmse$sb, col = "orange")
lines(100*expo, A.rmse$nl, col = "green")
lines(100*expo, A.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# WP
plot(100*expo, WP.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "WP RMSE (m)", ylab = "WP.rmse (m)",
     ylim = c(0,4))
lines(100*expo, WP.rmse$sb, col = "orange")
lines(100*expo, WP.rmse$nl, col = "green")
lines(100*expo, WP.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A0
plot(100*expo, A0.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)",
     main = "A0.rmse (m2)", ylab = "A0.rmse (m)", ylim = c(0,2000))
lines(100*expo, A0.rmse$sb, col = "orange")
lines(100*expo, A0.rmse$nl, col = "green")
lines(100*expo, A0.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# ------------------------------------------------------------------------------------------------
# Make plots of median z0, A0 profiles with error bounds

par(mfrow = c(2,2))

k <- 16 # exposure level
l.pred <- data.frame(r = 1:nr, 
                     z0.med = z0.l.med[,k], z0.lower = z0.l.lower[,k], z0.upper = z0.l.upper[,k],
                     A0.med = A0.l.med[,k], A0.lower = A0.l.lower[,k], A0.upper = A0.l.upper[,k],
                     s0.med = s0.l.med[,k], s0.lower = s0.l.lower[,k], s0.upper = s0.l.upper[,k],
                     WP.med = WP.l.med[,k], WP.lower = WP.l.lower[,k], WP.upper = WP.l.upper[,k],
                     A.med = A.l.med[,k], A.lower = A.l.lower[,k], A.upper = A.l.upper[,k])
sb.pred <- data.frame(r = 1:nr, 
                      z0.med = z0.sb.med[,k], z0.lower = z0.sb.lower[,k], z0.upper = z0.sb.upper[,k],
                      A0.med = A0.sb.med[,k], A0.lower = A0.sb.lower[,k], A0.upper = A0.sb.upper[,k],
                      s0.med = s0.sb.med[,k], s0.lower = s0.sb.lower[,k], s0.upper = s0.sb.upper[,k],
                      WP.med = WP.sb.med[,k], WP.lower = WP.sb.lower[,k], WP.upper = WP.sb.upper[,k],
                      A.med = A.sb.med[,k], A.lower = A.sb.lower[,k], A.upper = A.sb.upper[,k])
nl.pred <- data.frame(r = 1:nr, 
                      z0.med = z0.nl.med[,k], z0.lower = z0.nl.lower[,k], z0.upper = z0.nl.upper[,k],
                      A0.med = A0.nl.med[,k], A0.lower = A0.nl.lower[,k], A0.upper = A0.nl.upper[,k],
                      s0.med = s0.nl.med[,k], s0.lower = s0.nl.lower[,k], s0.upper = s0.nl.upper[,k],
                      WP.med = WP.nl.med[,k], WP.lower = WP.nl.lower[,k], WP.upper = WP.nl.upper[,k],
                      A.med = A.nl.med[,k], A.lower = A.nl.lower[,k], A.upper = A.nl.upper[,k])
nlsb.pred <- data.frame(r = 1:nr, 
                        z0.med = z0.nlsb.med[,k], z0.lower = z0.nlsb.lower[,k], z0.upper = z0.nlsb.upper[,k],
                        A0.med = A0.nlsb.med[,k], A0.lower = A0.nlsb.lower[,k], A0.upper = A0.nlsb.upper[,k],
                        s0.med = s0.nlsb.med[,k], s0.lower = s0.nlsb.lower[,k], s0.upper = s0.nlsb.upper[,k],
                        WP.med = WP.nlsb.med[,k], WP.lower = WP.nlsb.lower[,k], WP.upper = WP.nlsb.upper[,k],
                        A.med = A.nlsb.med[,k], A.lower = A.nlsb.lower[,k], A.upper = A.nlsb.upper[,k])

# z0
plot(z0.true.ra, main = paste("Median z0 (m) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(120,140), lwd = 1.5,
     xlab = "reach-average cross section", ylab = "z0")
lines(nl.pred$z0.lower, col = "lightgreen", lty = 2)
lines(nl.pred$z0.med, col = "green")
lines(nl.pred$z0.upper, col = "lightgreen", lty = 2)
lines(nlsb.pred$z0.lower, col = "lightblue", lty = 2)
lines(nlsb.pred$z0.med, col = "blue")
lines(nlsb.pred$z0.upper, col = "lightblue", lty = 2)
legend("bottomright", legend = c("True","NL","NLSB"),
       col = c("black","green","blue"), lwd = c(1,1,1))

# A0
plot(A0.true.ra[,k], main = paste("Median A0 at", expo[k]*100,"% exposure (sq. m)"), 
     type = "l", ylim = c(0,500), lwd = 1.5,
     xlab = "Cross section", ylab = "A0")
lines(nl.pred$A0.lower, col = "lightgreen", lty = 2)
lines(nl.pred$A0.med, col = "green")
lines(nl.pred$A0.upper, col = "lightgreen", lty = 2)
lines(nlsb.pred$A0.lower, col = "lightblue", lty = 2)
lines(nlsb.pred$A0.med, col = "blue")
lines(nlsb.pred$A0.upper, col = "lightblue", lty = 2)
legend("bottomright", legend = c("True","NL","NLSB"),
       col = c("black", "green","blue"), lwd = c(1,1,1))

# -----------------------------------------------------------------------------------------------------------------------------------
# Calculate average spread for A0 and z0 predictions

k <- 16
# mean(z0.l.upper[,k] - z0.l.lower[,k], na.rm = TRUE)
# mean(z0.sb.upper[,k] - z0.sb.lower[,k], na.rm = TRUE)
# # mean(z0.sbm.upper[,k] - z0.sbm.lower[,k], na.rm = TRUE)
# mean(z0.nl.upper[,k] - z0.nl.lower[,k], na.rm = TRUE)
# mean(z0.nlsb.upper[,k] - z0.nlsb.lower[,k], na.rm = TRUE)

# z0 IQR
z0.l.iqr <- apply(z0.l, c(1,2), IQR, na.rm = TRUE)
z0.sb.iqr <- apply(z0.sb, c(1,2), IQR, na.rm = TRUE)
z0.nl.iqr <- apply(z0.nl, c(1,2), IQR, na.rm = TRUE)
z0.nlsb.iqr <- apply(z0.nlsb, c(1,2), IQR, na.rm = TRUE)

# A0 IQR
A0.l.iqr <- apply(A0.l, c(1,2), IQR, na.rm = TRUE)
A0.sb.iqr <- apply(A0.sb, c(1,2), IQR, na.rm = TRUE)
A0.nl.iqr <- apply(A0.nl, c(1,2), IQR, na.rm = TRUE)
A0.nlsb.iqr <- apply(A0.nlsb, c(1,2), IQR, na.rm = TRUE)

# z0 sd
z0.l.sd <- apply(z0.l, c(1,2), sd, na.rm = TRUE)
z0.sb.sd <- apply(z0.sb, c(1,2), sd, na.rm = TRUE)
z0.nl.sd <- apply(z0.nl, c(1,2), sd, na.rm = TRUE)
z0.nlsb.sd <- apply(z0.nlsb, c(1,2), sd, na.rm = TRUE)

# A0 sd
A0.l.sd <- apply(A0.l, c(1,2), sd, na.rm = TRUE)
A0.sb.sd <- apply(A0.sb, c(1,2), sd, na.rm = TRUE)
A0.nl.sd <- apply(A0.nl, c(1,2), sd, na.rm = TRUE)
A0.nlsb.sd <- apply(A0.nlsb, c(1,2), sd, na.rm = TRUE)

# calculate mean IQR and sd for several different exposure levels
meanIQR <- vector(length = 4, "list") # first do for z0, then for A0
meanSD <- vector(length = 4, "list")
names(meanIQR) <- c("l","sb","nl","nlsb")
names(meanSD) <- c("l","sb","nl","nlsb")
for (k in 1:n_exp_levels)
{
  meanIQR$l[k] <- mean(A0.l.iqr[,k], na.rm = TRUE)
  meanIQR$sb[k] <- mean(A0.sb.iqr[,k], na.rm = TRUE)
  meanIQR$nl[k] <- mean(A0.nl.iqr[,k], na.rm = TRUE)
  meanIQR$nlsb[k] <- mean(A0.nlsb.iqr[,k], na.rm = TRUE)
  
  meanSD$l[k] <- mean(A0.l.sd[,k], na.rm = TRUE)
  meanSD$sb[k] <- mean(A0.sb.sd[,k], na.rm = TRUE)
  meanSD$nl[k] <- mean(A0.nl.sd[,k], na.rm = TRUE)
  meanSD$nlsb[k] <- mean(A0.nlsb.sd[,k], na.rm = TRUE)
}

# Check for normality: the IQR should be 1.34898*SD for a normally distributed RV
326*1.34898 # many of the predicted hydraulic variables appear Gaussian, based on this test.
# Should also do visual inspection

hist()
hist(A0.l[,,8])
hist(A0.nlsb[,16,])
IQR(A0.nlsb[,16,], na.rm = TRUE)
sd(A0.nlsb[,16,], na.rm = TRUE)

# -----------------------------------------------------------------------------------------------------------------------------------
# Plot prediction error histograms


# tmp <- 
# hist(tmp[,k,], main = paste("z0 error (m) at", expo[k]*100,"% exposure"), xlab = "z0 prediction error")
# sd(tmp[,k,])
# IQR(tmp[,k,])

# z0 prediction error
par(mfrow = c(4,4))
k <- 16
xl <- -4
xu <- 2
hist((z0.l - z0.true.ra)[,k,], 
     main = paste("z0 error (m) at", expo[k]*100,"% exposure, linear method"), 
     xlab = "z0 prediction error", "fd", xlim = c(xl,xu))
hist((z0.sb - z0.true.ra)[,k,], 
     main = paste("z0 error (m) at", expo[k]*100,"% exposure, SB method"), 
     xlab = "z0 prediction error", "fd", xlim = c(xl,xu))
hist((z0.nl - z0.true.ra)[,k,], 
     main = paste("z0 error (m) at", expo[k]*100,"% exposure, NL method"), 
     xlab = "z0 prediction error", "fd", xlim = c(xl,xu))
hist((z0.nlsb - z0.true.ra)[,k,], 
     main = paste("z0 error (m) at", expo[k]*100,"% exposure, NLSB method"), 
     xlab = "z0 prediction error", "fd", xlim = c(xl,xu))

# A0 prediction error
par(mfrow = c(4,4))
k <- 16
xl <- -100
xu <- 400
hist((A0.l - A0.true.ra[,k])[,k,], 
     main = paste("A0 error (sq. m) at", expo[k]*100,"% exposure, linear method"), 
     xlab = "A0 prediction error", "fd", xlim = c(xl,xu))
hist((A0.sb - A0.true.ra[,k])[,k,], 
     main = paste("A0 error (sq. m) at", expo[k]*100,"% exposure, SB method"), 
     xlab = "A0 prediction error", "fd", xlim = c(xl,xu))
hist((A0.nl - A0.true.ra[,k])[,k,], 
     main = paste("A0 error (sq. m) at", expo[k]*100,"% exposure, NL method"), 
     xlab = "A0 prediction error", "fd", xlim = c(xl,xu))
hist((A0.nlsb - A0.true.ra[,k])[,k,], 
     main = paste("A0 error (sq. m) at", expo[k]*100,"% exposure, NLSB method"), 
     xlab = "A0 prediction error", "fd", xlim = c(xl,xu))

# ------------------------------------------------------------------------------------------

# Figure N
#' Plot predicted vs. actual
#' 
#' Generates predicted vs. actual figure
#' Takes the median of the ensemble as the prediction

pred <- c(median(z0.l[1,k,]), median(z0.l[2,k,]), median(z0.l[3,k,]))
pred <- c(median(z0.sb[1,k,]), median(z0.sb[2,k,]), median(z0.sb[3,k,]))
pred <- c(median(z0.nl[1,k,]), median(z0.nl[2,k,]), median(z0.nl[3,k,]))
pred <- c(median(z0.nlsb[1,k,]), median(z0.nlsb[2,k,]), median(z0.nlsb[3,k,]))

plot(z0.true.ra, pred, 
     xlab = "true", ylab = "pred", asp = 1, 
     xlim = c(133,137), ylim = c(133,137), 
     main = "z0, 95% exposure, nlsb")
abline(0,1)

# Slope comparisons

k <- 19
summary(s0.l[2 , k, ])
hist(s0.l[2 , k, ])

r <- 1
s0.true.ra[r]
median(s0.l[r , k, ])
summary(s0.sb[r , k, ])
summary(s0.nl[r , k, ])
hist(s0.nlsb[r , k, ])

# ------------------------------------------------------------------------------------------------
# SCRAP 

# it would be good to take off outlier cross sections - say, those with bankfull parameters that are outliers
ra.xs <- vector(length = nr, "list")
# make this efficient later
sumw <- 0
sumWSE <- 0 
for (r in 1:2000)
{
  sumw <- sumw + xWSEw[[r]]$w
  sumWSE <- sumWSE + xWSEw[[r]]$WSE
}
ra.xs[[1]] <- data.frame(WSE = sumWSE/2000, w = sumw/2000)

sumw <- 0
sumWSE <- 0 
for (r in 2001:4000)
{
  sumw <- sumw + xWSEw[[r]]$w
  sumWSE <- sumWSE + xWSEw[[r]]$WSE
}
ra.xs[[2]] <- data.frame(WSE = sumWSE/2000, w = sumw/2000)

sumw <- 0
sumWSE <- 0 
for (r in 4001:5773)
{
  sumw <- sumw + xWSEw[[r]]$w
  sumWSE <- sumWSE + xWSEw[[r]]$WSE
}
# 5773-4001+1
ra.xs[[3]] <- data.frame(WSE = sumWSE/1773, w = sumw/1773)

plot(WSE~w, ra.xs[[1]])
plot(WSE~w, ra.xs[[2]])
plot(WSE~w, ra.xs[[3]])

###----

# --------------------------------------------------------------------------------------------------------------------------------------------
# Plot for presentation (one cross section)
r <- 1
WSEw_obs1 <- observe(rWSEw[[r]], exposure = 0.5)
# WSEw_obs1 <- observe(rWSEw[[r]], exposure = 0.5, sd_w = 0, sd_wse = 0)
plot(WSE~w, WSEw_obs1, xlim = c(0, 600), ylim = c(131,148), 
     main = paste("Reach averaged cross section", r, "at 50% exposure"),
     xlab = "w (m)", ylab = "WSE (m)")
lines(WSE~w, rWSEw[[r]])
# points(0, rWSEw[[1]]$WSE[1], col = "black", pch = 8)

# Fit models
lf1 <- fit_linear(WSEw_obs1)
sb1 <- fit_slopebreak(WSEw_obs1, multiple_breaks = FALSE, continuity = TRUE)
nl1 <- fit_nonlinear(WSEw_obs1)
nlsb1 <- fit_nlsb(WSEw_obs1)

# plot linear
lines(WSEw_obs1$w, predict(lf1), col = "red")
w.vals <- seq(0, min(WSEw_obs1$w), 1)
lines(w.vals, predict(lf1, newdata = data.frame(w=w.vals)), 
      col = "red", lty = 2)
# points(0, predict(lf1, newdata = data.frame(w=0)), col = "red", pch = 19)

# plot slope break
nn <- length(WSEw_obs1$w)
sb1.ind <- attributes(sb1)$sb.ind
lines(WSEw_obs1$w[1:sb1.ind], predict(sb1[[1]]), col = "orange")
lines(WSEw_obs1$w[(sb1.ind):nn], predict(sb1[[2]]), col = "orange")
lines(w.vals, predict(sb1[[1]], newdata = data.frame(w=w.vals)), 
      col = "orange", lty = 2)
# points(0, predict(sb1[[1]], newdata = data.frame(w=0)), col = "orange", pch = 19)

# plot nl
lines(WSEw_obs1$w, predict(nl1), col = "green")
lines(w.vals, predict(nl1, newdata = data.frame(w=w.vals)), 
      col = "green", lty = 2)
# points(0, predict(nl1, newdata = data.frame(w=0)), col = "green", pch = 19)

# plot nlsb
nlsb1.ind <- attributes(nlsb1)$sb.ind
lines(WSEw_obs1$w[1:nlsb1.ind], predict(nlsb1[[1]]), col = "blue")
lines(WSEw_obs1$w[(nlsb1.ind):nn], predict(nlsb1[[2]]), col = "blue")
lines(w.vals, predict(nlsb1[[1]], newdata = data.frame(w=w.vals)), 
      col = "blue", lty = 2)
# points(0, predict(nlsb1[[1]], newdata = data.frame(w=0)), col = "blue", pch = 19)

legend("topleft", legend = c("Data", "Linear","SB","NL","NLSB"), 
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# Predict hydraulic parameters for this one test case
A0.l <- array(dim = c(n_exp_levels, M))
A0.sb <- array(dim = c(n_exp_levels, M))
A0.nl <- array(dim = c(n_exp_levels, M))
A0.nlsb <- array(dim = c(n_exp_levels, M))
for (k in 10)
{
  for (m in 1:M)
  {
    w1 <- WSEw_obs1$w[1] # take the first value, supposing it is the minimum value
    h1 <- WSEw_obs1$WSE[1]
    A0.l[k,m] <- calc_model_A0(lf1, type = "linear")
    A0.sb[k,m] <- calc_model_A0(sb1, type = "sb", pos.only = FALSE)
    A0.nl[k,m] <- calc_model_A0(nl1, type = "nl", w1 = w1, h1 = h1, pos.only = FALSE)
    A0.nlsb[k,m] <- calc_model_A0(nlsb1, type = "nlsb", w1 = w1, h1 = h1, pos.only = FALSE)
  }
}

# Try to produce a negative value:
M <- 100
A0.nl <- array(dim = c(n_exp_levels, M))
bflag <-0
for (k in 1:n_exp_levels)
{
  for (m in 1:M)
  {
    WSEw_obs1 <- observe(rWSEw[[r]], exposure = expo[k])
    nl1 <- fit_nonlinear(WSEw_obs1)
    w1 <- WSEw_obs1$w[1] # take the first value, supposing it is the minimum value
    h1 <- WSEw_obs1$WSE[1]
    A0.nl[k,m] <- calc_model_A0(nl1, type = "nl", w1 = w1, h1 = h1, pos.only = FALSE)
    if (A0.nl[k,m] < 0)
    {
      bflag <- 1
      break
    }
  }
  if (bflag==1)
  {
    break
  }
}

# [1,24]
# [1,66]

# WSEw_obs1 <- observe(rWSEw[[r]], exposure = expo[1])


# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

#' # How come there are negative A0 predictions?
#' 
#' load(file.path(exp_dir, "pred_lf_bu.rda"))
#' length(pred_lf)
#' pred_lf[[1]]$A0[k,] # find a negative area
#' 
#' summary(A0.l)
#' which(A0.l<0, arr.ind = TRUE)
#' 
#' A0.l[1,1,1] # r, k, m
#' A0.l[1,3,156]
#' 
#' dim(A0.l)
#' 
#' summary(A0.l[1,1,])
#' 
#' r <- 1
#' k <- 1
#' A0.true.ra[r,k]
#' 
#' pred_lf[[1]]$A0[1,1]
#' 
#' WSEw_obs <- observe(rWSEw[[1]], exposure = 0.05, sd_w = 0, sd_wse = 0)
#' WSEw_obs <- observe(rWSEw[[1]], exposure = 0.05)
#' plot(WSE~w, WSEw_obs)
#' lf1 <- fit_linear(WSEw_obs)
#' calc_model_A0(lf1, type = "linear")
#' # abline(lf1)
#' 
#' #' If the slope is negative, then the A0 prediction will be negative. In some cases, especially when the error is relatively large, 
#' #' the model might describe a relationship where WSE decreases with width, which is physically inconsistent. 
#' #' Instead, the fitted model should be constrained to be positive only, or thrown out and redone with a different random error.

load(file.path(exp_dir, "pred_nlsb_bu.rda"))
length(pred_lf)
pred_nl[[1]]$A0[16,] # find a negative area

r <- 1

A0.true.ra[r,k]

k <- 16
WSEw_obs <- observe(rWSEw[[3]], exposure = k*0.05, sd_w = 0, sd_wse = 0)
# WSEw_obs <- observe(rWSEw[[3]], exposure = 1, sd_w = 0, sd_wse = 0)
WSEw_obs <- observe(rWSEw[[3]], exposure = k*0.05)
nlsb1 <- fit_nlsb(WSEw_obs)
nlsb1

plot(WSEw_obs$WSE[1:10], predict(nlsb1[[1]]), asp = 1, xlab = "true", ylab = "pred")
abline(0,1)

plot(WSE~w, WSEw_obs)
lines(WSEw_obs$w, predict(nl1), col = "red")

calc_model_A0(nlsb1, w1 = w0.ra[1,k], h1 = h1.ra[1,k], type = "nlsb", pos.only = TRUE)

MM <- 3
A0 <- vector(length = MM)
for (mm in 1:MM)
{
  WSEw_obs <- observe(rWSEw[[3]], exposure = k*0.05)
  # print(WSEw_obs)
  if (class(nlsb1[[1]]) == "nls")
  {
    # nlsb1 <- fit_nlsb(WSEw_obs)
    # print(coef(nlsb1[[1]])[1])
    nlsb_name <- paste0("nlsb/nlsb_", "r_", mm, ".rds")
    nlsb1 <- readRDS(file.path(exp_dir, nl_name))
    A0[mm] <- calc_model_A0(nlsb1, w1 = w0.ra[1,k], h1 = h1.ra[1,k], type = "nlsb", pos.only = TRUE)
    print(mm)
  }
}


WSEw_obs
nl1
w1 <- min(WSEw_obs$w)
# h1 <- min(fitted(nl1)) # this chooses the model-predicted WSE, rather than observed WSE, so it is not the right measure.
h1 <- min(WSEw_obs$WSE)
nl_area(w1, h1, coef(nl1)[1], coef(nl1)[2], coef(nl1)[3])

t1 <- (137.3407 - coef(nl1)[1])*w1
t2 <- (coef(nl1)[2]/(coef(nl1)[3]+1))*w1^(coef(nl1)[3]+1)
area0 <- t1 - t2 # hopefully t2<t1...

# calc_model_A0(lf1, type = "linear")

# --------------------------------------------------------------------------------------------------------
# Determining why A0 is sometimes negative

r <- 1
nlsb_name <- paste0("nlsb/nlsb_", "r_", mm, ".rds")
nlsb1 <- readRDS(file.path(exp_dir, nl_name))
nk <- 10
nm <- 10
A0 <- vector(length = nk*nm)
ind <- 1
for (k in 1:nk)
{
  for (m in 1:nm)
  {
    A0[ind] <- calc_model_A0(model = nlsb1[[k]][[m]], w1 = w0.ra[1,k], h1 = h1.ra[1,k], type = "nl")
    ind <- ind + 1
  }
}

h1 <- h1.ra[1,k]
w1 = w0.ra[1,k]
t1 <- (h1 - coef(nlsb1[[k]][[m]])[1])*w1
t2 <- (coef(nlsb1[[k]][[m]])[2]/(coef(nlsb1[[k]][[m]])[3]+1))*w1^(coef(nlsb1[[k]][[m]])[3]+1)
area0 <- t1 - t2

# This model has a negative A0:

nlsb_name <- paste0("nlsb/nlsb_", "r_", r, ".rds")
nlsb1 <- readRDS(file.path(exp_dir, nlsb_name))
k <- 10
m <- 10
nlsb1 <- nlsb1[[k]][[m]]

# Here are the data used to fit this model
obsname <- paste0("obs/WSEw_obs_r_", r, ".rds")
WSEw_obs1 <- readRDS(file.path(exp_dir, obsname))
WSEw_obs1 <- WSEw_obs1[[k]][[m]]

w1 <- min(WSEw_obs1$w)
h1 <- min(WSEw_obs1$WSE)

calc_model_A0(model = nlsb1, w1 = w1, h1 = h1, type = "nlsb")

# Plot the fitted model

plot(WSE~w, WSEw_obs1, xlim = c(0, 600), ylim = c(131,144), 
     main = paste("Reach averaged cross section", r, "at 50% exposure"),
     xlab = "w (m)", ylab = "WSE (m)")
lines(WSE~w, rWSEw[[r]])

nn <- length(WSEw_obs1$w)
nlsb1.ind <- attributes(nlsb1)$sb.ind
lines(WSEw_obs1$w[1:nlsb1.ind], predict(nlsb1[[1]]), col = "blue")
lines(WSEw_obs1$w[(nlsb1.ind):nn], predict(nlsb1[[2]]), col = "blue")
lines(w.vals, predict(nlsb1[[1]], newdata = data.frame(w.vals= seq(0, min(WSEw_obs1$w), 1))), 
      col = "blue", lty = 2)



