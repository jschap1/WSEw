# Analyze cross section characteristics
#
# Use this to generate tables and figures for presentation

setwd("/Users/jschap/Documents/Research/SWOTBATH")
savename = "p26_channel_params.rda"
load(savename)

summary(wbf)
summary(dbf.max)
summary(dbf)
summary(A)
summary(s)

sd(wbf)
sd(dbf.max)
sd(dbf)
sd(A)
sd(s)

n.xs <- length(wbf)

# This is all performed by characterize_channel(), so it is duplicate here.
# par(mfrow = c(4,1))
# plot(dist_downstream/1000, dbf, type = "l", 
#      main = "maximum bankfull depth (m)", xlab = "distance downstream (km)",
#      ylab = "dbf_max")
# plot(dist_downstream/1000, dbf.max, type = "l", 
#      main = "maximum bankfull depth (m)", xlab = "distance downstream (km)",
#      ylab = "dbf_max")
# plot(dist_downstream/1000, wbf, type = "l", 
#      main = "bankfull width (m)", xlab = "distance downstream (km)",
#      ylab = "wbf")
# plot(dist_downstream/1000, b.min, type = "l", 
#      main = "minimum bed elevation (m)", xlab = "distance downstream (km)",
#      ylab = "b_min")
# plot(dist_downstream/1000, A, type = "l", 
#      main = "bankfull flow area (sq. m)", 
#      xlab = "distance downstream (km)",
#      ylab = "A")
# 
# par(mfrow = c(2,1))
# plot(dist_downstream/1000, s, type = "l", 
#      main = "shape parameter", 
#      xlab = "distance downstream (km)",
#      ylab = "s")
# abline(1,0, col = "gray")
# abline(2,0, col = "gray")
# hist(s, main = "shape parameters", "fd", col = "lightgray")
  

