# Multiple slope break method (batch)
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
numlocations <- vector(length = length(exposure))

for (i in 1:length(exposure))
{
  
  # Initialize vectors
  sb.location <- vector(length = nseg)
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
    WSEw <- na.omit(WSEw)
    
    # Test if optimal location for multiple slope break method
    win <- 4
    if (dim(WSEw)[1] <= 2*win){next} # error check
    
    # Test for first slope break
    sb.ind <- test_slope_break(WSEw$WSE, WSEw$w, thres = 0.015, window = win, m = TRUE)
    
    sb.add <- 1
    # Test for additional slope breaks
    if (!is.null(sb.ind))
    {
      while(!is.null(sb.add))
      {
        if (length(WSEw$WSE[1:sb.ind[length(sb.ind)]])<=2*win) 
        {
          break
        }
        sb.add <- test_slope_break(WSEw$WSE[1:sb.ind[length(sb.ind)]], 
                              WSEw$w[1:sb.ind[length(sb.ind)]], 
                              thres = 0.015, window = win, m = TRUE)
        sb.ind <- append(sb.ind, sb.add)
      }
      if (sb.ind == 1)
      {
        next
      }
    }
    
    if (!is.null(sb.ind))
    {
      
      # If location is optimal, estimate WSE0
      sb.location[seg] <- TRUE
      lf <- lm(WSE~w, data = WSEw[1:sb.ind[length(sb.ind)],])
      summary(lf)
      r2 <- cor(WSEw)[1,2]^2
      WSE0[seg] <- lf$coefficients[1]
      
      # Calculate error of the estimate (absolute error w.r.t. b.min)
      WSE0_error[seg] <- WSE0[seg] - b.min[seg]
      
    }
  }
  numlocations[i] <- sum(sb.location)
  RMSE[i] <- sqrt(sum(WSE0_error^2))
}

par(mfrow=c(1,1))
hist(WSE0_error[sb.location])
summary(WSE0_error[sb.location])

# ------------------------------------------------------------------------
# Slope break method summary

par(mfrow = c(1,1), mar = c(5,5,2,5))
df <- data.frame(nl = numlocations, RMSE = RMSE, exp = exposure)

with(df, plot(100*exp, RMSE, main = "Depth RMSE for Multiple Slope Break Method",
              xlab = "Channel exposure (%)", 
              ylab = "RMSE (m)", pch = 19, lwd = 2, col = "red"))

par(new = TRUE)
with(df, plot(100*exp, nl, pch = 25, col = "blue", lwd = 2, 
              axes = FALSE, xlab = NA, ylab = NA))

axis(side = 4)
mtext(side = 4, line = 3, "Optimal locations")
legend("top", col = c("red", "blue"), pch = c(19, 25),
       legend=c("RMSE","#locations"))

# ------------------------------------------------------------------------
# Plotting 

# Cross section
par(mfrow=c(1,1))
plot(cross.section$x, cross.section$b, type = "l")

# Cross section and WSE-w relationship
par(mfrow = c(1,2))
plot(cross.section$x, cross.section$b, xlab = "x (m)", ylab = "WSE (m)", 
     main = paste("Pool 21, Segment", seg))
lines(cross.section$x, f1(tdist[[seg]]))
plot(w, WSE, xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship")

# ------------------------------------------------------------------------
# Etc.

plot_WSEw <- function(seg)
  # Plots WSE-w relationship
{
  cross.section <- data.frame(x = tdist[[seg]], d = d[[seg]], b = b[[seg]]) 
  f1 <- approxfun(cross.section$x, y = cross.section$b, method="linear",
                  yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)
  WSE <- seq(b.min[seg], b.min[seg] + dbf[seg], by = interval*dbf[seg])
  w <- get_width2(WSE, f1, x.max = max(cross.section$x), delx = 1)
  par(mfrow = c(1,1))
  plot(w, WSE, xlab = "width (m)", ylab = "WSE (m)", 
       main = "Width-WSE relationship")
  return(NULL)
}

for (seg in 1:100)
{
  png(paste0("WSEw", seg, ".png"))
  plot_WSEw(seg)
  dev.off()
}

# seg = 13, 25 are good test cases for multiple slope break method

# For final slope break, test whether the trend is linear
res <- check_below_sb(WSEw$WSE, WSEw$w, 
                      sb.ind = sb.ind[length(sb.ind)], thres = 0.015)
if (is.null(res))
{
  sb.ind <- NULL
}

