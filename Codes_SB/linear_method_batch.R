# Linear method (batch)
# Following Mersel et al. (2013)
#
# Inputs: 
# transects.depth (depth)
# b (bed elevation)

rm(list=ls())
setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Slope_Break")
load("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Transects/pool_21_5m.rda")

library(rootSolve)
source("sb_functions.R")
opar <- par()

interval <- 0.05 # interval for entering in with WSE values
exposure <- seq(0.1, 1, by = 0.1) # channel exposure

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

# Initialize vectors
RMSE <- vector(length = length(exposure))
numlinearlocations <- vector(length = length(exposure))

for (i in 1:length(exposure))
{
  
  # Initialize vectors
  linear.location <- vector(length = nseg)
  WSE0 <- vector(length = nseg)
  WSE0_error <- vector(length = nseg)
  
for (seg in 1:nseg)
{
  
  # Estimate WSE-w relationship
  cross.section <- data.frame(x = tdist[[seg]], d = d[[seg]], b = b[[seg]]) 
  
  # Fit a linear spline to the observed channel geometry
  f1 <- approxfun(cross.section$x, y = cross.section$b, method="linear",
                  yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)
  
  # Make a vector of WSE(t) values to enter in with
  
  if (dbf[seg] < 0){next} # error check
  
  WSE <- seq(b.min[seg], b.min[seg] + dbf[seg], by = interval*dbf[seg])
  
  # Estimate the flow width for each WSE value
  w <- get_width2(WSE, f1, x.max = max(cross.section$x), delx = 1)
  
  # Put WSE-w info into a data frame
  WSEw <- data.frame(WSE = WSE, w = w)
  
  # Remove unobserved channel geometry
  if (exposure[i]<1)
  {
    unobserved.ind <- which(WSEw$w/wbf[seg] < (1-exposure[i]))
    WSEw <- WSEw[-unobserved.ind,]
  }
  
  if (dim(WSEw)[1] <= 1){next} # error check
  if (all(is.na(WSEw$w))) {next} # error check
  
  # Test if optimal location for linear method
  s1 <- test_linear(WSEw$WSE, WSEw$w, thres = 0.015)
  if (!is.null(s1))
  {
    linear.location[seg] <- TRUE
    # If linear method is optimal, use the linear method
    lf <- lm(WSE~w, data = WSEw)
    summary(lf)
    r2 <- cor(WSEw)[1,2]^2
    WSE0[seg] <- lf$coefficients[1]
    
    # Calculate error of the estimate (absolute error w.r.t. b.min)
    WSE0_error[seg] <- WSE0[seg] - b.min[seg]
    
  }
}
  numlinearlocations[i] <- sum(linear.location)
  RMSE[i] <- sqrt(sum(WSE0_error^2))
}


hist(WSE0_error[linear.location])
summary(WSE0_error[linear.location])

# ------------------------------------------------------------------------
# Linear method summary

par(mfrow = c(1,1), mar = c(5,5,2,5))
df <- data.frame(nl = numlinearlocations, RMSE = RMSE, exp = exposure)

with(df, plot(100*exp, RMSE, main = "Channel depth RMSE for Linear Method",
              xlab = "Channel exposure (%)", log = "y", 
              ylab = "RMSE (m)", pch = 19, lwd = 2, col = "red"))

par(new = TRUE)
with(df, plot(100*exp, nl, pch = 25, col = "blue", lwd = 2, 
              axes = FALSE, xlab = NA, ylab = NA))

axis(side = 4)
mtext(side = 4, line = 3, "Optimal locations")
legend("topright", col = c("red", "blue"), pch = c(19, 25),
        legend=c("RMSE","#locations"))
     
