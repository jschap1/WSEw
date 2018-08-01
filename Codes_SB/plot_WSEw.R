# Plot all WSEw relationships and cross section shapes

rm(list=ls())
setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Slope_Break/Codes_07272018")
load("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections/Slope_Break/Data/processed_data_p21_500m.rda")
opar <- par()

nseg <- length(cross.sections)
for (seg in 1:nseg)
{
  png(paste0("p21_seg_", seg, ".png"))
  par(mfrow = c(1,2))
  plot(cross.sections[[seg]]$x, cross.sections[[seg]]$b, xlab = "x (m)", ylab = "WSE (m)", 
       main = paste("Pool 21, Segment", seg), type = "l", ylim = c(b.min[seg], b.min[seg] + dbf[seg]))
  plot(xWSEw[[seg]]$w, xWSEw[[seg]]$WSE, xlab = "width (m)", ylab = "WSE (m)", 
       main = "Width-WSE relationship")
  dev.off()
}

# -------------------------------------------------------------------------------------
# Plot all WSEw relationships on one plot
par(mfrow = c(1,1))
seg <- 1
plot(xWSEw[[seg]]$w, xWSEw[[seg]]$WSE, xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship", type = "l", lwd = 0.5, xlim = c(0, 1500), ylim = c(130,145))
for (seg in 2:nseg)
{
  lines(xWSEw[[seg]]$w, xWSEw[[seg]]$WSE, lwd = 0.5)
}

par(mfrow = c(1,1))
seg <- 1
plot(rWSEw[[seg]]$w, rWSEw[[seg]]$WSE, xlab = "width (m)", ylab = "WSE (m)", 
     main = "Width-WSE relationship", type = "l", lwd = 0.5, xlim = c(0, 1000), ylim = c(134,143))
for (seg in 2:nseg)
{
  lines(rWSEw[[seg]]$w, rWSEw[[seg]]$WSE, lwd = 0.5)
}





# -------------------------------------------------------------------------------------
# Finding the correlation length for a cross section

cross.sections[[seg]]$x
cross.sections[[seg]]$b

df <- data.frame(x = cross.sections[[seg]]$x,b = cross.sections[[seg]]$b)

library(gstat)
# g <- gstat(id = "b", formula = b, locations = ~x, data = df)

# Make 2D data frame
channel.pix <- vector(length=nseg)
for (seg in 1:nseg)
{
  channel.pix[seg] <- length(cross.sections[[seg]]$x)
}

a <- array(dim = c(max(channel.pix), nseg, 3))

grd <- expand.grid(1:291, 1:nseg)
names(grd) <- c("x","seg")

a <- vector(length = dim(grd)[1])
ind <- 1
for (seg in 1:nseg)
{
  for (x in 1:291)
  {
    if (!is.na(cross.sections[[seg]]$b[x]))
    {
      a[ind] <- cross.sections[[seg]]$b[x]
    } else
    {
      a[ind] <- NA
    }
    ind <- ind + 1
  }
}

b <- cbind(grd, a)
names(b) <- c("x","s","b")

library(geoR)
g <- as.geodata(b)
par(mfrow = c(1,1), mar = c(5,5,2,5))
points(g, xlab = "cross distance", ylab = "distance downstream")

var1 <- variog(g)
plot(var1)
var1 <- variog(g, trend = "1st")
plot(var1)
var1 <- variog(g, estimator.type = "modulus")
plot(var1)