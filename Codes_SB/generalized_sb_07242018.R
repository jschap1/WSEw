# Generalized linear/slope break method
# Extending the slope break method of Mersel et al. (2013) to allow nonlinear fits
# In development 7/24/2018
# Designed to run for a single segment
# A batch method will be developed as well, but it will be a separate script.
#
# Inputs: 
# transects.depth (depth)
# b (bed elevation)
# seg (segment number)

rm(list=ls())
setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Slope_Break")
load("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Transects/pool_21_5m.rda")
opar <- par()

library(rootSolve)
source("sb_functions.R")

# ------------------------------------------------------------------------
# Pre-processing

# Remove empty channels
rm.ind <- which(channel.pix<=10)
transects.depth <- transects.depth[-rm.ind]
b <- b[-rm.ind]
tdist <- tdist[-rm.ind]

# Pre-processing
nseg <- length(b)
d <- transects.depth
dbf <- unlist(lapply(d, max, na.rm=TRUE)) # bankfull depth (m)
wbf <- unlist(lapply(tdist, max, na.rm=TRUE)) # bankfull width (m)
b.min <- unlist(lapply(b, min))





w <- tdist
wbf <- unlist(lapply(tdist, max, na.rm=TRUE)) # bankfull width (m)

b
bbf <- lapply(b, max, na.rm=TRUE)

d <- transects.depth
dbf <- unlist(lapply(d, max, na.rm=TRUE)) # bankfull depth (m)

seg <- 70


cross.section <- preprocess(seg, w, wbf, b.smooth, bbf, rescale = TRUE)
names(cross.section) <- c("norm.x","norm.d")
cross.section <- na.omit(cross.section)
par(opar)
plot(cross.section, type = "l")

nseg <- length(b)

b.min <- unlist(lapply(b, min))

# ------------------------------------------------------------------------
# Estimate WSE-w relationship

# Fit a linear spline to the measured channel geometry
f1 <- approxfun(cross.section$norm.x, y = cross.section$norm.d, method="linear",
                yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)

# Make a vector of values to enter in with
d.in <- seq(0, 1, by = 0.01)

# Estimate the flow width for each WSE value
w <- get_width2(d.in, f1, x.max = 1, delx = 0.001)

norm.wd <- data.frame(norm.w = w, norm.d = d.in)

# Remove unobserved part of cross section
#unobserved.ind <- which(norm.wd$norm.w < 0.2)
#norm.wd <- norm.wd[-unobserved.ind,]

# Remove the vertical part of the cross section
norm.wd <- na.omit(norm.wd)

# Make plot
par(mfrow = c(1,2))
plot(cross.section$norm.x, cross.section$norm.d, xlab = "Normalized x", ylab = "Normalized depth", 
     main = paste("Pool 21, Segment", seg), type = "l")
plot(norm.d~norm.w, data = norm.wd, xlab = "Normalized width", ylab = "Normalized depth", 
     main = "Width-depth relationship", ylim = c(0,1))

# ------------------------------------------------------------------------
# Check if there are one or more shape breaks

# Test for first shape break
window <- 30
sb.thres <- 0.3
sb.ind <- test_gen_sb(d.in, w, sb.thres <- sb.thres, window = window)

# Test for a second shape break
sb.add <- test_gen_sb(d.in[1:sb.ind[length(sb.ind)]], w[1:sb.ind[length(sb.ind)]], sb.thres = sb.thres)

# Test for additional shape breaks
while(!is.null(sb.add))
{
  if (length(d.in[1:sb.ind[length(sb.ind)]])<=window) 
  {
    break
  }
  sb.add <- test_gen_sb(d.in[1:sb.ind[length(sb.ind)]], w[1:sb.ind[length(sb.ind)]], sb.thres = sb.thres)
  sb.ind <- append(sb.ind, sb.add)
}

# sb.ind <- c(9) # manually entering shape break locations
n.sb <- length(sb.ind) # number of slope breaks

# Plot shape break locations
par(opar)
plot(w,d.in, main = "shape break location")
points(w[sb.ind], d.in[sb.ind], col = "red", pch = 19, lwd = 5)

# Locations can sometimes be just a little off

# ------------------------------------------------------------------------
# Plot the normalized cross section fits

# Split the data according to shape breaks (automate later)
nn <- dim(norm.wd)[1]
wd <- vector(length = n.sb+1, "list")
wd[[1]] <- norm.wd[1:sb.ind[n.sb],]
wd[[2]] <- norm.wd[(sb.ind[n.sb]+1):(sb.ind[n.sb-1]-1),]
wd[[3]] <- norm.wd[(sb.ind[n.sb-(n.sb-1)]):nn,]

wd[[1]] <- norm.wd[1:3,]
wd[[2]] <- norm.wd[4:11,]
wd[[3]] <- norm.wd[12:21,]

wd[[1]] <- norm.wd[1:7,]
wd[[2]] <- norm.wd[7:8,]
wd[[3]] <- norm.wd[9:18,]

# Fit curves to each set of points
s <- vector(length = n.sb+1)
fit.nls <- vector(length = n.sb+1, "list")
norm.w.pred <- vector(length = n.sb+1, "list")
for (k in 1:(n.sb+1))
{
  fit.nls[[k]] <- nls(norm.w~norm.d^(1/s),start = list(s = 1), data=wd[[k]])
  s[k] <- coef(fit.nls[[k]])
  norm.w.pred[[k]] <- predict(fit.nls[[k]])
  lines(norm.w.pred[[k]], wd[[k]]$norm.d, col = "red")
}

# ------------------------------------------------------------------------
# Use the fitted WSE-w relationship to estimate unobserved bottom depth

d.est <- dbf*(w/wbf)^s[1] # not sure the shape approach is good...

# Might be simpler to stop thinking about it in terms of shape parameter, and just fit the cross-sections
# WSE-w with a general power model, without normalizing. 

# ------------------------------------------------------------------------
# Make a plot showing the cross section, the width-WSE relationship, and the fitted cross section

b.min <- min(b.smooth[[seg]], na.rm = TRUE)
d <- transects.depth
dbf <- unlist(lapply(d, max, na.rm=TRUE)) # bankfull depth (m)

w.pred <- vector(length = n.sb+1, "list")
WSE <- vector(length = n.sb+1, "list")
for (k in 1:(n.sb+1))
{
  w.pred[[k]] <- wbf[[seg]]*norm.w.pred[[k]] # predicted width
  WSE[[k]] <- b.min + wd[[k]]$norm.d*dbf[seg] # WSE inputs
}

w.unnorm <- wbf[[seg]]*w # unnormalized flow width

cross.section.unnorm <- data.frame(x = tdist[[seg]], d = d[[seg]], b = b.smooth[[seg]]) # geometry of cross-section

par(mfrow = c(1,2))
plot(cross.section.unnorm$x, cross.section.unnorm$b, type = "l",
     main = "Pool 21, Segment 8", ylab = "WSE (m)", xlab = "x (m)")

##########
# Add plot of idealized cross section
# There are extra bumps in the plot bc I have not forced continuity at the breakpoints

x.mid <- mean(cross.section.unnorm$x)
x.pred.left <- x.mid - unlist(w.pred)/2
x.pred.right <- x.mid + unlist(w.pred)/2
b.pred <- unlist(WSE)

lines(na.omit(x.pred.left), b.pred, col = "red")
lines(na.omit(x.pred.right), b.pred, col = "red")

##########

plot(na.omit(w.unnorm), unlist(WSE), xlim = c(0, wbf[[seg]]), 
     ylim = c(b.min, max(cross.section.unnorm$b, na.rm = TRUE)), main = "WSE-width relationship",
     xlab = "width (m)", ylab = "WSE (m)")

lines(w.pred[[1]], WSE[[1]], col = "red")
lines(w.pred[[2]], WSE[[2]], col = "red")
lines(w.pred[[3]], WSE[[3]], col = "red")

# ------------------------------------------------------------------------
# Compare area and wetted perimeter of the original and fitted cross sections (if time)


# ------------------------------------------------------------------------
# SCRAP

# plot(dat$width, dat$depth, xlab = "x (m)", ylab = "WSE (m)", 
main = paste("Pool 21, Segment", seg), ylim = c(min(WSE),max(WSE)))
# lines(dat$x, f1(tdist[[seg]]))
plot(w*wbf[[seg]], WSE, xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship")
lines(w.pred, WSE[1:sb.ind], xlab = "width (m)", ylab = "WSE (m)", col = "cyan")
legend("topleft", legend = paste("s.low =", round(s.low, digits = 3)))











# fit a curve to the uppermost WSE-w pairs
wd.high <- data.frame(w = w[(sb.ind+1):nn], d = d.in[(sb.ind+1):nn])
fit.nls.high <- nls(w~d^(1/s),start = list(s = 1), data=wd.high)
s.high <- coef(fit.nls)

# fit a curve to the next set of WSE-w pairs
dat.low <- data.frame(w = w[1:sb.ind], d = d.in[1:sb.ind])
fit.nls.low <- nls(w~d^(1/s),start = list(s = 1), data=dat.low)
s.low <- coef(fit.nls.low)

# Plot results
par(mfrow = c(1,2))
plot(dat$width, dat$depth, xlab = "Normalized x", ylab = "Normalized depth", 
     main = paste("Pool 21, Segment", seg), type = "l")
plot(w, d.in, xlab = "Normalized width", ylab = "Normalized depth", 
     main = "Width-depth relationship")

lines(predict(fit.nls.high), dat.high$d, col = "cyan")
lines(predict(fit.nls.low), dat.low$d, col = "red")


# ------------------------------------------------------------------------
# Estimate the bottom depth using the lowermost fitted curve



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




# ------------------------------------------------------------------------
# Generalized "linear" method

# Fit s value
datfd <- data.frame(w = w, d = d.in)
fit.nls <- nls(w~d^(1/s),start = list(s = 1), data=datfd)
s <- coef(fit.nls)

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

