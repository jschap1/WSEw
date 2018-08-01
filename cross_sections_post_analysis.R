# Plot histogram of shape parameters from all pools
#
# June 28, 2018

rm(list=ls())
setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections")
opar <- par()

pools <- c(4,21) # pools analyzed

s.all <- vector()
# Load shape parameters
for (p in pools)
{
  # load(paste0("Transects/s_p", p, ".rda")) # loads s
  load(paste0("Transects/s_p", p, "_w2_nls.rda")) # loads s
  s.all <- c(s.all,s)
}
s.all <- s.all[-which(is.na(s.all))]

par(mfrow=c(1,1))
hist(s.all, breaks = seq(0,ceiling(max(s.all)),length.out = 10), col = "darkblue",
     main = "Shape parameters for UMRB Pools 4, 21", xlab = "s")
summary(s.all)

