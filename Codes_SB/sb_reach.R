# 7/25/2018
# Reach averaged slope break method

rm(list=ls())
setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Slope_Break")
load("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Transects/pool_21_500m.rda")
source("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Slope_Break/Codes_07272018/sb_functions.R")
opar <- par()

result <- calc_WSEw(interval = 0.05, dx = 1, 
                    geometry_data = "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Transects/pool_21_500m.rda") # should really choose a number that corresponds with SWOT measurements in a year
xWSEw <- result$WSEw
cross.sections <- result$xs

rWSEw <- reach_avg(xWSEw, l = 10000, res = 500)

exposure <- seq(0.1, 1, by = 0.1) # channel exposure

# Remove empty channels
if (any(channel.pix<=10))
{
  rm.ind <- which(channel.pix<=10)
  transects.depth <- transects.depth[-rm.ind]
  b <- b[-rm.ind]
  tdist <- tdist[-rm.ind]
}

# Pre-processing
nseg <- length(b)
d <- transects.depth
dbf <- unlist(lapply(d, max, na.rm=TRUE)) # bankfull depth (m)
b.min <- unlist(lapply(b, min))

# returns max width, as long as none of the reaches have larger WSE than max(width)
wbf <- unlist(lapply(xWSEw, max))

save(xWSEw, rWSEw, cross.sections, b.min, wbf, dbf, file = "processed_data_p21_500m.rda")

# b.min is the average minimum bed elevation for the reach
nr <- length(rWSEw)
b.min.xs <- unlist(lapply(b, min, na.rm=T)) # min bed elev for a cross section
b.min <- vector(length = nr) # reach average min bed elev 
l <- 10000
res <- 5
pix <- l/res
for (r in 1:nr)
{
  b.min[r] <- mean(b.min.xs[r:(r+(pix-1))])
}

# Initialize vectors
RMSE <- vector(length = length(exposure))
numlocations <- vector(length = length(exposure))

for (i in 1:length(exposure))
{
  
  # Initialize vectors
  sb.location <- vector(length = nr)
  WSE0 <- vector(length = nr)
  WSE0_error <- vector(length = nr)
  
  for (r in 1:nr)
  {
    
    WSEw <- rWSEw[[r]] # WSE-w relationship for a reach
    nn <- dim(WSEw)[1]
    
    # Remove unobserved channel geometry
    if (exposure[i]<1)
    {
      unobserved.ind <- which(WSEw$w/wbf[r] < (1-exposure[i]))
      WSEw <- WSEw[-unobserved.ind,]
    }
    
    # Test if optimal location for slope break method
    win <- 4
    if (dim(WSEw)[1] <= 2*win){next} # error check
    sb.ind <- test_slope_break(WSEw$WSE, WSEw$w, thres = 0.015, window = win)
    if (!is.null(sb.ind))
    {
      # If location is optimal, estimate WSE0
      sb.location[r] <- TRUE
      lf.upper <- lm(WSE~w, data = WSEw[sb.ind:nn,])
      lf <- lm(WSE~w, data = WSEw[1:sb.ind,])
      summary(lf)
      r2 <- cor(WSEw)[1,2]^2
      WSE0[r] <- lf$coefficients[1]
      
      #z0 <- lf$coefficients[1]
      #s <- lf$coefficients[2]
      
      # WSE0 estimate with standard error
      # Reverse the roles of w and WSE bc WSE variable is deterministic here
      # And the standard error can be interpreted this way
      #sb.location[seg] <- TRUE
      #lf <- lm(w~WSE, data = WSEw[1:sb.ind,])
      #z0 <- -lf$coefficients[1]/lf$coefficients[2]
      #s <- 1/lf$coefficients[2]
      # First order approximation of the variance, assuming independent coefficient estimates
      #var_z0 <- ((36268^2)/(262.7^2))*((753^2)/(36268^2)+(5.415^2)/(262.7^2))
      #sd_z0 <- sqrt(var_z0)
      
      # Calculate error of the estimate (absolute error w.r.t. b.min)
      WSE0_error[r] <- WSE0[r] - b.min[r]
      
    }
  }
  numlocations[i] <- sum(sb.location)
  RMSE[i] <- sqrt(sum(WSE0_error^2))
}

hist(WSE0_error[sb.location])
summary(WSE0_error[sb.location])

# ------------------------------------------------------------------------
# Slope break method summary

par(mfrow = c(1,1), mar = c(5,5,2,5))
df <- data.frame(nl = numlocations, RMSE = RMSE, exp = exposure)

with(df, plot(100*exp, RMSE, main = "Channel depth RMSE for Slope Break Method",
              xlab = "Channel exposure (%)", 
              ylab = "RMSE (m)", pch = 19, lwd = 2, col = "red"))

par(new = TRUE)
with(df, plot(100*exp, nl, pch = 25, col = "blue", lwd = 2, 
              axes = FALSE, xlab = NA, ylab = NA))

axis(side = 4)
mtext(side = 4, line = 3, "Optimal locations")
legend("topright", col = c("red", "blue"), pch = c(19, 25),
       legend=c("RMSE","#locations"))

# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# Calculate standard error of estimate (manually)

# Assumes deterministic predictor variables (w)
N <- dim(WSEw[1:sb.ind,])[1]
X <- cbind(rep(1,N), WSEw[1:sb.ind,]$w)
SSE <- sum(lf$residuals^2)
C <- (1/(N-2))*SSE*solve(t(X)%*%X) # variance-covariance matrix
SE1 <- sqrt(C[1,1])
SE2 <- sqrt(C[2,2]) # Compare to summary(lf)
# Note that there is non-negligible cross-covariance
# So there is a relationship between slope and intercept
# When slope is steep, intercept is small, and vice versa

# ------------------------------------------------------------------------
# Make a plot showing the width-WSE relationship, and the idealized cross section

r <- 4031

par(mfrow = c(2,1))
plot(WSE~w, data = rWSEw[[r]], xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship")
predWSE <- predict(lf, newdata = lf$model)
predWSE.upper <- predict(lf.upper, newdata = lf.upper$model)
lines(lf$model$w, predWSE, col="blue")
lines(lf.upper$model$w, predWSE.upper, col="blue")

plot(lf$model$w/2, predWSE, xlim = c(-wbf[r], wbf[r]), type = "l",
     main = "Idealized cross section", col = "red", ylim = c(b.min[r], max(predWSE.upper)))
lines(-lf$model$w/2, predWSE, col = "red")
lines(lf.upper$model$w/2, predWSE.upper, col = "red")
lines(-lf.upper$model$w/2, predWSE.upper, col = "red")

# ------------------------------------------------------------------------
# Nonlinear fit

r <- 4031

par(opar)
plot(WSE~w, data = rWSEw[[r]], xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship")

z0.true <- 134.44

# How does accuracy change with exposure?
exposure <- seq(0.1, 1, by = 0.1)

wbf <- max(rWSEw[[r]]$w)

# Remove unobserved channel geometry
if (exposure[i]<1)
{
  unobserved.ind <- which(rWSEw[[r]]$w/wbf < (1-exposure[i]))
  rWSEw[[r]] <- rWSEw[[r]][-unobserved.ind,]
}

fit.nls <- nls(WSE~a+b*w^(c), start = list(a=2, b=2, c=2), data=rWSEw[[r]])
wse.pred <- predict(fit.nls)
lines(rWSEw[[r]]$w, wse.pred,  col = "red")

coefs <- coefficients(fit.nls)
z0 <- coefs[1]
err <- z0 - z0.true

# ------------------------------------------------------------------------
# Final for abstract - developing the method

r <- 4031

rWSEw.copy <- rWSEw # make a backup copy
rWSEw  <- rWSEw.copy

exposure <- 0.5 # half the banks are exposed
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

# Monte Carlo approach to SWOT measurement uncertainty

M <- 1e3 # number of MC simulations
set.seed(704753262)

nn <- dim(rWSEw[[r]])[1] # number of observations

#e_WSE <- rnorm(nn, mean = 0, sd = 0.25) # standard deviation of 0.25 m
#e_w <- rnorm(nn, mean = 1, sd = 0.25) # cv = 0.25

par(mfrow = c(2,1))
hist(e_WSE, breaks = 40, main = "WSE error, 1e3 samples", col = "darkblue")
hist(e_w, breaks = 40, main = "width error, 1e3 samples", col = "darkblue")
summary(e_w)

par(opar)
plot(WSE~w, data = rWSEw.copy[[r]], xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship", lty = 1, type = "l", lwd = 3, ylim = c(130, 146), xlim = c(0, 1000))

z0 <- vector(length = M)
sb.ind.debug <- vector(length = M)
for (m in 1:M)
{
  # corrupt observations
  
  e_WSE <- rnorm(nn, mean = 0, sd = 0.25) # standard deviation of 0.25 m
  e_w <- rnorm(nn, mean = 1, sd = 0.25) # cv = 0.25
  
  WSEw <- rWSEw[[r]]
  WSEw$WSE <- rWSEw[[r]]$WSE + e_WSE 
  WSEw$w <- (e_w)*rWSEw[[r]]$w
  
  # points(WSE~w, data = WSEw, cex = 0.3, col = "cyan")
  
  # find slope break
  sb.ind <- test_slope_break(WSEw$WSE, WSEw$w, thres = 0.5, window = 10, m = FALSE)
  if (is.null(sb.ind)) 
  {
    sb.ind.debug[m] <- NA
    next
  } else
  {
    sb.ind.debug[m] <- sb.ind
  }
  
  WSEw1 <- WSEw[1:sb.ind,] # points below the slope break
  lf1 <- lm(WSE~w, data = WSEw1)
  z0[m] <- lf1$coefficients[1]
  
  lines(seq(0, 800, by = 1), predict(lf1, newdata = data.frame(w = seq(0, 800, by = 1))), col = "darkgreen", lwd = 0.5)
  
}

z0[z0==0] <- NA
hist(z0, breaks = 40)

# calculate some statistics
z0.true <- rWSEw.copy[[r]]$WSE[1]
mse <- (1/(M-sum(is.na(z0))))*sum((z0-rep(z0.true, M))^2, na.rm = TRUE)
bias <- mean(z0, na.rm = TRUE) - z0.true
variance <- mse-bias
stddev <- variance^0.5

plot(WSE~w, data = rWSEw.copy[[r]], xlab = "width (m)", ylab = "WSE (m)", 
      main = "Width-WSE relationship", lty = 1, type = "l", lwd = 3)

lines(WSE~w, data = rWSEw[[r]], xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship", lty = 1, type = "l", lwd = 3)

lines(WSEw1$w, predict(lf1), col = "red")

hist(z0, breaks = 40)
summary(sb.ind)





