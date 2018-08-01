# More on fitting cross sections
# Goals:
# Identify "nice" cross sections (check)
# Fit all the nice cross sections (check)

rm(list=ls())
setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Cross_Sections")
load("Transects/pool_21_w2.rda")
opar <- par()

library(rootSolve)

source("fitting_sections_functions.R")
source("fit_shape.R")

# ------------------------------------------------------------------------
# Pre-processing

nseg <- length(transects.depth)
w <- tdist
wbf <- widths
b <- b.smooth
bbf <- lapply(b, max, na.rm=TRUE)
d <- d.smooth # depth
dbf <- lapply(d, max, na.rm=TRUE)

seg <- 8
dat <- preprocess(seg, w, wbf, b, bbf, rescale = TRUE)
dat <- na.omit(dat)
plot(dat)

# ------------------------------------------------------------------------
# Compute symmetry SSE, a statistic I am using to measure symmetry

nseg <- length(b)
SSE <- vector(length = nseg)
for (seg in 1:nseg)
{
  if (channel.pix[seg] <= 1) # zero-width channels
  { # error handling
    SSE[seg] <- NA
    next
  } else if (max(b[[seg]], na.rm=TRUE) == min(b[[seg]], na.rm=TRUE)) # flat cross sections
  {
    SSE[seg] <- NA
    next    
  }
  dat <- preprocess(seg, w, wbf, b, bbf, rescale = TRUE)
  dat <- na.omit(dat)
  SSE[seg] <- calc_symm_SSE(dat)
}

# ------------------------------------------------------------------------
# Identify symmetric cross sections

# Visual inspection method
# plot all transects so I can see which are symmetric
for (seg in 1:nseg)
{
  if (channel.pix[seg] == 0) # zero-width channels
  { # error handling
    next
  } else if (max(b[[seg]], na.rm=TRUE) == min(b[[seg]], na.rm=TRUE)) # flat cross sections
  {
    next    
  }
  dat <- preprocess(seg, w, wbf, b, bbf, rescale = TRUE)
  png(paste0("Figures/Pool5/Transects/transect", seg, "w2.png"))
  plot(dat, main = paste("Cross Section", seg, "symmetry SSE =", round(SSE[seg],2)))
  dev.off()
}

# Use the lowest quantile as a threshold
summary(SSE)
hist(SSE, breaks = seq(0,ceiling((max(SSE, na.rm=T))), length.out = 15))
thres <- as.numeric(summary(SSE)[2])
symm.seg <- which(SSE<=thres) # numbers of the segments that are relatively symmetric

# ------------------------------------------------------------------------
# Identify single-channel sections

# Finds the local minima in the channel cross section
# Uses smoothing techniques described at https://rpubs.com/mengxu/peak_detection

nseg <- length(b)
nchannels <- vector(length = nseg)
for (seg in 1:nseg)
{
  if (channel.pix[seg] <= 1) # zero-width channels
  { # error handling
    next
  } else if (max(b[[seg]], na.rm=TRUE) == min(b[[seg]], na.rm=TRUE)) # flat cross sections
  {
    next    
  }  
  dat <- preprocess(seg, w, wbf, b, bbf, rescale = TRUE)
  dat <- na.omit(dat)
  x <- dat$width
  y <- -dat$depth
  peaks <- test(w = 5, span = 0.5) # tuning parameters
  nchannels[seg] <- length(peaks$x)
}

chan1.seg <- which(nchannels==1)

# ------------------------------------------------------------------------
# Return the "best" sections for fitting

nice.ind <- symm.seg[symm.seg %in% chan1.seg]

# ------------------------------------------------------------------------
# Get shape parameters for these cross sections

s <- vector(length = nseg)
for (seg in nice.ind)
{
  dat <- preprocess(seg, w, wbf, b, bbf, rescale = TRUE)
  dat <- na.omit(dat)
  # s[seg] <- fit_shape(dat)
  s[seg] <- fit_nls(dat)
}
s[-nice.ind] <- NA

save(s, file = "Transects/s_p5_w2_nls.rda")

hist(s, main = "shape parameters for fitted cross sections")
length(nice.ind)
summary(s)

# Plot them
for (seg in nice.ind)
{
  dat <- preprocess(seg, w, wbf, b, bbf, rescale = TRUE)
  png(paste0("Figures/Pool4/transects_nls/w2_", seg, ".png")) 
  plot(dat, main = paste0("Segment", seg, ", s = ", round(s[seg],2) ))
  # also add the modeled cross section to the plot
  dev.off()
}

# It would be great to mark the locations of these cross sections on the map of pool 21 (use auto_transects)
# Also would be good to show diagnostics for the fits

# -------------------------------------------------------------------- 
# Modeled cross section - make a plot of the true and modeled cross section

seg <- 8

# indexing issues mean this only works for even lengths...
c.len <- channel.pix[seg]
  
if ((c.len)%%2 == 0)
{
  wl <- w[[seg]][1:(c.len/2)]
  wr <- w[[seg]][(c.len/2+1):c.len]
} else
{
  wl <- w[[seg]][1:floor(c.len/2)]
  wr <- w[[seg]][(ceiling(c.len/2)):c.len]
}

br <- min(b[[seg]], na.rm=T) + bbf[[seg]]*(wl/wbf[[seg]])^s[seg] #[seg]
bl <- rev(br)

if ((c.len)%%2 == 0)
{
  b.model <- c(bl, br)
} else
{
  b.model <- c(bl, NA, br)
}

par(mfrow=c(1,1))
png(paste0("Figures/Pool5/Fitted_Channels/channel_", seg, "_w2.png"))
plot(c(wl, wr), b[[seg]], xlab = "x", ylab = "bed elevation (m)", 
     main = paste("Segment", seg, "s=", round(s[seg],2)), lty = 1, type = "l")
lines(c(wl, wr), b.model, type = "l", col = "blue", lty = 2)
legend("top", col = c("black","blue"), lty = c(1,2), legend = c("Actual","Fitted"))
dev.off()


