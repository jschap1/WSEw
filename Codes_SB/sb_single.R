# Linear and/or slope break method
# Following Mersel et al. (2013)
#
# Inputs: 
# transects.depth (depth)
# b (bed elevation)
# seg (segment number)

rm(list=ls())
# setwd("C:/Users/Jacob/Desktop/cross_sections_R_codes_07192018")
setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Slope_Break")
load("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Transects/pool_21_5m.rda")
opar <- par()

library(rootSolve)
source("sb_functions.R")

#source("fitting_sections_functions.R")
#source("fit_shape.R")

# ------------------------------------------------------------------------
# Pre-processing

d <- transects.depth
seg <- 5783
cross.section <- data.frame(x = tdist[[seg]], d = d[[seg]], b = b[[seg]]) # geometry of cross-section
opar <- par()
par(mfrow=c(1,1))
plot(cross.section$x, cross.section$b)

# ------------------------------------------------------------------------
# Estimate WSE-w relationship

# Fit a linear spline to the observed channel geometry
f1 <- approxfun(cross.section$x, y = cross.section$b, method="linear",
                yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)

# Make a vector of WSE(t) values to enter in with
dbf <- unlist(lapply(d, max, na.rm=TRUE)) # bankfull depth (m)
b.min <- min(b[[seg]])
WSE <- seq(b.min, b.min + dbf[[seg]], by = 0.05*dbf[[seg]])

# Estimate the flow width for each WSE value
w <- get_width2(WSE, f1, x.max = max(cross.section$x))

# Make plot
par(mfrow = c(1,2))
plot(cross.section$x, cross.section$b, xlab = "x (m)", ylab = "WSE (m)", 
     main = paste("Pool 21, Segment", seg), ylim = c(min(WSE),max(WSE)))
lines(cross.section$x, f1(tdist[[seg]]))
plot(w, WSE, xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship")

# ------------------------------------------------------------------------
# Test if the linear method will work

s1 <- test_linear(WSE,w, thres=0.015) # returns first difference slope estimates or NULL
if (!is.null(s1))
{
  plot(w, WSE, main = "WSE-w")
  plot(w, s1, main = "slope")
}

# ------------------------------------------------------------------------
# Estimate the bottom depth using the linear method

WSEw <- data.frame(WSE = WSE, w = w)
lf <- lm(WSE~w, data = WSEw)

summary(lf)
r2 <- cor(WSEw)[1,2]^2

par(mfrow = c(1,1))
plot(w, WSE, main = "Linear")
lines(WSEw$w, predict(lf, newdata = WSEw))
legend("topleft", legend = paste("r^2 =", round(r2, digits = 4)))

# ------------------------------------------------------------------------
# Test if the slope break method will work

sb.ind <- test_slope_break(WSE,w)
if (!is.null(sb.ind))
{
  plot(w, WSE, main = "WSE-w")
  points(w[sb.ind], WSE[sb.ind], col = "red", pch = 19, lwd = 5)
}

# ------------------------------------------------------------------------
# Estimate the bottom depth using the slope break method

# fit a line to the lowermost WSE-w pairs
dat.low <- data.frame(WSE = WSE[1:sb.ind], w = w[1:sb.ind])
lf.low <- lm(WSE~w, data = dat.low)

summary(lf.low)
r2 <- cor(dat.low)[1,2]^2

par(mfrow = c(1,1))
plot(w, WSE, main = "Slope Break")
lines(dat.low$w, predict(lf.low, newdata = dat.low))
points(w[sb.ind], WSE[sb.ind], col = "red", pch = 2, lwd = 3)
legend("topleft", legend = paste("r^2 =", round(r2, digits = 4)))

