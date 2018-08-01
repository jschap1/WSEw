# Generalized slope break method
# Adaptation of Mersel et al. (2013) for a general exponential cross section shape
# In development 7/16/2018
#
# Normalizing widths and depths allows us to use the shape parameter as an unknown
# Using pre-specified shape parameter means we can use unnormalized widths and depths to more faithfully follow Mersel's approach

rm(list=ls())
# setwd("C:/Users/Jacob/Desktop/cross_sections_R_codes_07192018")
setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections")
load("Transects/pool_21_w2.rda")
opar <- par()

library(rootSolve)
source("gen_slope_break_functions.R")

#source("fitting_sections_functions.R")
#source("fit_shape.R")



# ------------------------------------------------------------------------
# Pre-processing

d <- transects.depth; b
seg <- 35
dat <- data.frame(x = tdist[[seg]], d = d[[seg]], b = b[[seg]])
par(mfrow=c(1,1))
plot(dat$x, dat$b)

# Fit a linear spline to the measured channel geometry
f1 <- approxfun(dat$x, y = dat$b, method="linear",
                yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)

# Make a vector of WSE(t) values to enter in with
dbf <- unlist(lapply(transects.depth, max, na.rm=TRUE)) # bankfull depth (m)
b.min <- min(b[[seg]]) # b[[seg]][bf.ind[seg]]
WSE <- seq(b.min, b.min + dbf[[seg]], by = 0.05*dbf[[seg]])

# Estimate the flow width for each WSE value
w <- get_width2(WSE, f1, x.max = max(dat$x))

# Make plot
opar <- par()
par(mfrow = c(1,2))
plot(dat$x, dat$b, xlab = "x (m)", ylab = "WSE (m)", 
     main = paste("Pool 21, Segment", seg), ylim = c(min(WSE),max(WSE)))
lines(dat$x, f1(tdist[[seg]]))
plot(w, WSE, xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship")

# ------------------------------------------------------------------------
# Test if the linear method will work (eventually, loop over transects)

nn <- length(WSE)
s1 <- vector(length = nn)
for (k in 1:(nn-1))
{
  s1[k] <- (WSE[k+1]-WSE[k])/(w[k+1]-w[k]) # forward difference
}
s1[nn] <- s1[n-1]
# remove infinite values (where width did not change)
s1[s1==Inf] <- NA

plot(w, WSE, main = "WSE-w")
plot(w, s1, main = "slope")

thres <- 0.015 # determine through trial and error 
maxdiff <- max(s1, na.rm = TRUE) - min(s1, na.rm=TRUE)

if (maxdiff < thres)
{
  print("Use the linear method")
} else
{
  print("Do not use the linear method")
}

# ------------------------------------------------------------------------
# Test if the slope break method will work

# calculate mean of the slope for the highest few water surface
window <- 4
s1.high <- s1[(nn-(window-1)):nn] 
s.bar <- mean(s1.high, na.rm = TRUE)

# compare to s1 at lower water surfaces to find a slope break
for (j in 1:(window+1))
{
  if (s1[nn-(window-1)-j] < 0.3*s.bar)
  {
    sb.ind <- nn-(window-1)-j
    print("There is a slope break here")
    break
  } else
  {
    s.bar <- mean(s.bar, s1[nn-(window-1)-j], na.rm = TRUE)
  }
}

plot(w,WSE, main = "slope break location")
points(w[sb.ind], WSE[sb.ind], col = "red")

# Check that behavior is linear below the slope break

s2 <- s1[1:sb.ind]
maxdiff <- max(s2, na.rm = TRUE) - min(s1, na.rm=TRUE)

if (maxdiff < thres)
{
  print("Use the slope break method")
} else
{
  print("Do not use the slope break method")
}

# ------------------------------------------------------------------------
# Estimate the bottom depth using the slope break method

# fit a line to the lowermost WSE-w pairs
dat.low <- data.frame(WSE = WSE[1:sb.ind], w = w[1:sb.ind])
f2 <- lm(WSE~w, data = dat.low)

summary(f2)
r2 <- cor(dat.low)[1,2]^2

par(mfrow = c(1,1))
# plot(dat.low$w, dat.low$WSE, main = "lower WSE-w pairs")
plot(w, WSE, main = "Slope Break")
lines(dat.low$w, predict(f2, newdata = dat.low))
points(w[sb.ind], WSE[sb.ind], col = "red", pch = 2, lwd = 3)
legend("topleft", legend = paste("r^2 =", round(r2, digits = 4)))

# ------------------------------------------------------------------------
# Generalized "linear" method

w <- tdist
wbf <- widths
bbf <- lapply(b, max, na.rm=TRUE)

seg <- 35
dat <- preprocess(seg, w, wbf, b, bbf, rescale = TRUE)
dat <- na.omit(dat)
plot(dat)

# Fit a linear spline to the measured channel geometry
f1 <- approxfun(dat$width, y = dat$depth, method="linear",
                yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)

# Make a vector of values to enter in with
d.in <- seq(0, 1, by = 0.05)

# Estimate the flow width for each WSE value
w <- get_width2(d.in, f1, x.max = 1, delx = 0.01)

# Fit s value
datfd <- data.frame(w = w, d = d.in)
fit.nls <- nls(w~d^(1/s),start = list(s = 1), data=datfd)
s <- coef(fit.nls)

# Make plot
par(mfrow = c(1,2))
plot(dat$w, dat$d, xlab = "Normalized x", ylab = "Normalized depth", 
     main = paste("Pool 21, Segment", seg), ylim = c(min(WSE),max(WSE)), type = "l")
plot(d.in~w, data = datfd, xlab = "Normalized width", ylab = "Normalized depth", 
     main = "Width-depth relationship")
lines(predict(fit.nls), datfd$d, col = "cyan")

# Use fitted model to make plots in terms of WSE and w
b.min <- min(b[[seg]])
WSE <- b.min + d.in*dbf[seg]
w.pred <- wbf[[seg]]*predict(fit.nls)

par(mfrow = c(1,2))
plot(dat$x, dat$b, xlab = "x (m)", ylab = "WSE (m)", 
     main = paste("Pool 21, Segment", seg), ylim = c(min(WSE),max(WSE)))
lines(dat$x, f1(tdist[[seg]]))
plot(w*wbf[[seg]], WSE, xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship")
lines(w.pred, WSE, xlab = "width (m)", ylab = "WSE (m)", col = "cyan")
legend("topleft", legend = paste("s =", round(s, digits = 3)))


# ------------------------------------------------------------------------
# Check for "slope breaks"

# criterion: does s parameter change by a lot as you fit successive sets of points?
# If so, then there is a "slope break"
# More of a "shape break" really

window <- 4
high.ind <- (nn-3):nn
datfd <- data.frame(w = w[high.ind], d = d.in[high.ind])
fit.nls <- nls(w~d^(1/s),start = list(s = 1), data=datfd)
s.high <- coef(fit.nls)

# compare to s at lower water surfaces to find a "slope break"
for (j in 1:(nn-2*window))
{
  low.ind <- (nn-(2*window-1+j)):(nn-(window+j))
  datfd.low <- data.frame(w = w[low.ind], d = d.in[low.ind]) 
  fit.nls <- nls(w~d^(1/s),start = list(s = 1), data=datfd.low)
  s.low <- coef(fit.nls) 
  if (s.low < 0.3*s.high)
  {
    sb.ind <- low.ind[window]
    print("There is a slope break here")
    break
  } else
  {
    high.ind <- (nn-(window-1+j)):(nn)
    datfd <- data.frame(w = w[high.ind], d = d.in[high.ind])
    fit.nls <- nls(w~d^(1/s),start = list(s = 1), data=datfd)
    s.high <- coef(fit.nls)
  }
}

plot(w,d.in, main = "slope break location")
points(w[sb.ind], d.in[sb.ind], col = "red")

# Check that behavior is "linear" below the slope break (don't bother for now...)

#s2 <- s1[1:sb.ind]
#maxdiff <- max(s2, na.rm = TRUE) - min(s1, na.rm=TRUE)

if (maxdiff < thres)
{
  print("Use the slope break method")
} else
{
  print("Do not use the slope break method")
}

# ------------------------------------------------------------------------
# Do the nonlinear "slope break" method


# fit a curve to the uppermost WSE-w pairs
dat.high <- data.frame(w = w[(sb.ind+1):nn], d = d.in[(sb.ind+1):nn])
fit.nls.high <- nls(w~d^(1/s),start = list(s = 1), data=dat.high)
s.high <- coef(fit.nls)

# fit a curve to the lowermost WSE-w pairs
dat.low <- data.frame(w = w[1:sb.ind], d = d.in[1:sb.ind])
fit.nls.low <- nls(w~d^(1/s),start = list(s = 1), data=dat.low)
s.low <- coef(fit.nls)

# Plot results
par(mfrow = c(1,2))
plot(dat$width, dat$depth, xlab = "Normalized x", ylab = "Normalized depth", 
     main = paste("Pool 21, Segment", seg), type = "l")
plot(w, d.in, xlab = "Normalized width", ylab = "Normalized depth", 
     main = "Width-depth relationship")

lines(predict(fit.nls.high), dat.high$d, col = "cyan")
lines(predict(fit.nls.low), dat.low$d, col = "red")

# Use fitted model to make plots in terms of WSE and w
b.min <- min(b[[seg]])
WSE <- b.min + d.in*dbf[seg]
w.pred <- wbf[[seg]]*predict(fit.nls.low)

par(mfrow = c(1,2))
# plot(dat$width, dat$depth, xlab = "x (m)", ylab = "WSE (m)", 
     main = paste("Pool 21, Segment", seg), ylim = c(min(WSE),max(WSE)))
# lines(dat$x, f1(tdist[[seg]]))
plot(w*wbf[[seg]], WSE, xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship")
lines(w.pred, WSE[1:sb.ind], xlab = "width (m)", ylab = "WSE (m)", col = "cyan")
legend("topleft", legend = paste("s.low =", round(s.low, digits = 3)))


# ------------------------------------------------------------------------
# Do a fit for a cross section

seg <- 8
dat <- preprocess(seg, w, wbf, b, bbf, rescale = TRUE)
dat <- na.omit(dat)
plot(dat)

#fit_nls(dat) # fitting one function
#fit_nls2(dat) # fitting two functions


datfd <- get_w_WSE_relation(dat)

# Determine whether the linear or slope break method should be used
t <- test_method(datfd)

fit.nls <- nls(w~d^(1/s),start = list(s = 1),data=datfd)
s <- coef(fit.nls)

# Plot fit
par(mfrow=c(1,2))
plot(dat)
plot(d.in, flow_width, xlab = "normalized depth", ylab = "normalized width",
     main = paste("Shape parameter = ", round(s,2)), cex = 1, col = "cyan")
lines(datfd$d, predict(fit.nls))

# ------------------------------------------------------------------------
# For all cross sections
t <- vector(length = nseg)
for (seg in 1:nseg)
{
  dat <- preprocess(seg, w, wbf, b, bbf, rescale = TRUE)
  dat <- na.omit(dat)
  if (is.null(dim(dat))) {next}
  datfd <- get_w_WSE_relation(dat)
  t[seg] <- test_method(datfd)
}

seg <- 52 # or 35
dat <- preprocess(seg, w, wbf, b, bbf, rescale = TRUE)
dat <- na.omit(dat)
plot(dat)




