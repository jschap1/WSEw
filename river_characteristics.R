# Characterizing the river channel
# Adapt from old plotting functions to make some good plots 
# of h-w relationships and fits

setwd("/Users/jschap/Desktop/Cross_Sections")
opar <- par()
load("Data/Processed/processed_xs_data.rda")

library(WSEw)

n.xs <- length(xWSEw)

# Plot river characteristics
dist_downstream <- seq(5, n.xs*5, by = 5) # m

# bankfull width
wbf <- unlist(lapply(cross_sections$x, max))
plot(dist_downstream/1000, wbf, type = "l", 
     main = "bankfull width (m)", xlab = "distance downstream (km)",
     ylab = "wbf")
# plot(ra(wbf, 2000), type = "l")

# maximum bankfull depth
dbf.max <- unlist(lapply(cross_sections$d, max))
plot(dist_downstream/1000, dbf.max, type = "l", 
     main = "maximum bankfull depth (m)", xlab = "distance downstream (km)",
     ylab = "dbf_max")

# average bankfull depth
dbf <- unlist(lapply(cross_sections$d, mean))
plot(dist_downstream/1000, dbf, type = "l", 
     main = "maximum bankfull depth (m)", xlab = "distance downstream (km)",
     ylab = "dbf_max")

# minimum bed elevation
b.min <- unlist(lapply(cross_sections$b, min))
plot(dist_downstream/1000, b.min, type = "l", 
     main = "minimum bed elevation (m)", xlab = "distance downstream (km)",
     ylab = "b_min")

# shape parameter (using linearized fitting method to avoid sensitivity to initial guesses)
xdw <- vector(length = n.xs, "list")
A <- vector(length = n.xs)
for (r in 1:n.xs)
{
  width <- xWSEw[[r]]$w
  depth <- xWSEw[[r]]$WSE - b.min[r]
  A[r] <- calc_A_from_WSEw(xWSEw[[r]]) # bankfull flow area
  xdw[[r]] <- data.frame(w = width/wbf[r], d = depth/dbf.max[r]) # width-depth data
}

# perform fits
drop_zero <- TRUE
add_small <- FALSE
power.fits <- vector(length = n.xs, "list")
s <- vector(length = n.xs)
for (r in 1:n.xs)
{
  if (drop_zero) # need to modify the data to account for zero values
  {
    dw <- xdw[[r]][-1,]
  } else if (add_small) # this really messes up the fits, don't do this
  {
    dw <- xdw[[r]]
    dw$w[1] <- 1e-6
    dw$d[1] <- 1e-6
  }
  power.fits[[r]] <- lm(log(d) ~ log(w) + 0, data = dw)
  s[r] <- as.numeric(coef(power.fits[[r]]))
}

# power.fits[[r]]
hist(s)
summary(s)

# plot h-w for a cross section
plot(d~w, xdw[[r]])

# bankfull flow area
summary(A)
summary(ra(A,2000))
plot(A, type = "l")
plot(ra(A,2000), type = "l")
