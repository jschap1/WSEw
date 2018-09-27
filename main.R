# Executive file for WSEw project
# 
# Created 8/2/2018 JRS
# Revised 8/10/2018 JRS
# Revised 9/5/2018 JRS 
#   Reorganized for ease of use
#   Removed scrap code from end of document
#   Changed looping structure to avoid having to load all models into memory at once
# Revised 9/6/2018 JRS
#   Added Monte Carlo simulation capability. Can set M = 1 replicate for deterministic case.
#   Changed the way models are stored to avoid having to store too much data in RAM.
# Revised 9/13/2018
#   Implemented multicore processing with foreach package for model fitting
#   Removed timers
#   Put model-fitting functions in a separate R script
#   Updating plotting functions to use medians
#   Note: some parts of the code have not been updated, so be careful when running it
# Revised 9/17/2018
#   Implemented multicore processing for prediction
# Revised 9/21/2018
#   Removed SBM functionality, eliminated redundancies, vectorized loops
#   This is the optimization before running the code many times.
# Revised 9/26/2018
#   Changed reach averaging method, now uses 3, 10 km reaches for the 30 km study area (pool 21)
#   Started implementing SWOT data from gauge time series instead of even sampling method,
#   but currently still using even sampling.
#   Changed name of executive file to "main.R" without a "dN" on the end

# ------------------------------------------------------------------------------------------------

# Set up environment
library(foreach)
library(doMC)
library(raster)
library(strucchange)
library(minpack.lm)
library(WSEw)

setwd("/Users/jschap/Desktop/Cross_Sections")

library(devtools)
library(roxygen2)
update_WSEw()

# Experiment description
nr <- 3
reach_avg <- "10km"
spacing <- 5 # m
swot_sampling <- "even"
pool <- 21
err_type <- "no_err"
M <- 1 # number of replicates

# Make a directory to store results
exp_desc <- paste0("pool_", pool, "_ra_",reach_avg,"_nr_",nr,"_spacing_",spacing,"_sampling_",swot_sampling, "_", err_type, "_replicates_", M)
fits_dir <- "/Volumes/HD3/Cross_Sections" # directory for modeling outputs
exp_dir <- file.path(fits_dir, exp_desc) # directory for this experiment's outputs

if (!dir.exists(exp_dir))
{
  dir.create(exp_dir)
  dir.create(file.path(exp_dir, "lf"))
  dir.create(file.path(exp_dir, "sb"))
  dir.create(file.path(exp_dir, "nl"))
  dir.create(file.path(exp_dir, "nlsb"))
  dir.create(file.path(exp_dir, "obs"))
}

set.seed(704753262)
source("Codes/modelfit_par.R") # load the fitting functions
source("Codes/predict_par.R") # load the predicting functions

# ------------------------------------------------------------------------------------------------
# Load raw data

refWSE <- 470 # refWSE for pool 21 (ft)
refWSE <- refWSE/(39.37/12) #  convert to m

# Bathymetry
if (!file.exists("Data/p21_depth.tif"))
{
  bathy.dir <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC"
  bathy.name <- "bath_pool_21/bath_1999_p21"
  bathyfile <- file.path(bathy.dir, bathy.name)
  #levels(umesc)
  umesc <- raster(bathyfile)
  depth_5 <- as_numeric_raster(umesc, att = "DEPTH_M") # native resolution depths (5 m)
  writeRaster(depth_5, file = "Data/p21_depth.tif")
} else 
{
  depth_5 <- raster("Data/p21_depth.tif")
  # depth_50 <- aggregate(depth_5, fact = 10) # resample depth to 50 m resolution
}

# River centerline
riv.dir <- "/Users/jschap/Desktop/Cross_Sections/Data/Centerlines"
load(file.path(riv.dir, "centerline21.rda"))
riv <- centerline_p21

# Load USGS stream gauges in UMRB
gauges <- read.table("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/Stage/gauges_UMB_QC.txt")
names(gauges) <- c("id","lat","lon","V4")
coordinates(gauges) <- ~lon + lat # make SpatialPoints object
crs(gauges) <- "+init=epsg:4326"
gauges.utm <- spTransform(gauges, crs(riv))

# Plot the study area
plot(depth_5, main = "UMRB Pool 21", xlab = "easting", ylab = "northing")
lines(riv)
points(gauges.utm, pch = 19, col = "black", cex = 0.5)

# ------------------------------------------------------------------------------------------------
# Load processed cross section and WSE-w data

load("Data/Processed/processed_xs_data.rda")

# ------------------------------------------------------------------------------------------------
# Load fitted models (don't try to load them all at once)

lf <- readRDS(file.path(exp_dir, "lf.rds"))
sb <- readRDS(file.path(exp_dir, "sb.rds")) 

# ------------------------------------------------------------------------------------------------
# Load true hydraulic parameters

w0.ra <- readRDS(file.path(exp_dir, "w0_ra.rds"))
WP.true.xs <- readRDS(file.path(exp_dir, "WP_true_xs.rds"))
WP.true.ra <- readRDS(file.path(exp_dir, "WP_true_ra.rds")) 

A.true.xs <- readRDS(file.path(exp_dir, "A_true_xs.rds"))
A.true.ra <- readRDS(file.path(exp_dir, "A_true_ra.rds"))
A0.true.ra <- readRDS(file.path(exp_dir, "A0_true_ra.rds"))

s0.true.ra <- readRDS(file.path(exp_dir, "s0_true_ra.rds"))
s0.true.xs <- readRDS(file.path(exp_dir, "s0_true_xs.rds"))
z0.true.ra <- readRDS(file.path(exp_dir, "z0_true_ra.rds"))
z0.true.xs <- readRDS(file.path(exp_dir, "z0_true_xs.rds"))

# ------------------------------------------------------------------------------------------------
# Load observations
# WSEw_obs <- readRDS(file.path(exp_dir, "WSEw_obs.rds"))

# ------------------------------------------------------------------------------------------------
# Load predicted hydraulic parameters

load(file.path(exp_dir, "z0_pred.rda"))
load(file.path(exp_dir, "A_pred.rda"))
load(file.path(exp_dir, "WP_pred.rda"))
load(file.path(exp_dir, "A0_pred.rda"))
load(file.path(exp_dir, "s0_pred.rda"))

# ------------------------------------------------------------------------------------------------
# Load prediction statistics

load(file.path(exp_dir, "bias.rda"))
load(file.path(exp_dir, "rmse.rda"))

# ------------------------------------------------------------------------------------------------
# Compute cross section data from raw bathymetry

cross_sections <- auto_transects(section_length = 5, depth = depth_5, refWSE = refWSE, 
                                 savename = transects_name, makeplot = FALSE, riv = riv)
# (Takes about 4 hours at the highest possible resolution)

# Method 1: find corresponding flow width for WSE ranging from empty to bankfull conditions
xWSEw <- calc_WSEw(cross_sections, interval = 0.05, dx = 1) # number of data points depends on discretization

rWSEw <- reach_avg(xWSEw)
rWSEw_burr <- reach_avg(xWSEw)

save(cross_sections, xWSEw, rWSEw, file = file.path(exp_dir, "processed_xs_data_10km.rda"))

# Plot the observations
plot(WSE~w, rWSEw[[1]], main = "WSE-w sampling, three years")
points(WSE~w, rWSEw_burr[[1]], col="red", pch=19)
legend("topleft", legend = c("Even sampling","Burr sampling"), fill = c("black", "red"))

# ------------------------------------------------------------------------------------------------
# Make reach-average effective cross sections that are 10 km long each

# hist(unlist(lapply(cross_sections$x, length)))
xs.res <- resample_xs(cross_sections, n = 300)
xs.avg <- vector(length = 3, "list")
xs.avg[[1]] <- calc_mean_cross_section(xs.res[1:2000,,])
xs.avg[[2]] <- calc_mean_cross_section(xs.res[2001:4000,,])
xs.avg[[3]] <- calc_mean_cross_section(xs.res[4001:5773,,])

# ------------------------------------------------------------------------------------------------
# Calculate h-w relationship for the reach-average cross sections

# need to reformat
x <- vector("list", length = nr)
b <- vector("list", length = nr)
d <- vector("list", length = nr)
for (r in 1:nr)
{
  x[[r]] <- xs.avg[[r]]$x
  b[[r]] <- xs.avg[[r]]$b
  d[[r]] <- xs.avg[[r]]$d
}
cross_sections_avg <- list(x = x, b = b, d = d) 
rWSEw <- calc_WSEw(cross_sections_avg, interval = 0.05, dx = 1)

save(cross_sections, xWSEw, rWSEw, file = file.path(exp_dir, "processed_xs_data_10km.rda"))

par(mfrow = c(3,2))
plot(b~x, xs.avg[[1]], type = "l", 
     main = "Reach average cross section 1",
     xlim = c(0,800), ylim = c(135,144))
plot(WSE~w, rWSEw[[1]], xlab = "width (m)", ylab = "height (m)", 
     ylim = c(135, 144), main = "Height-width relationship", type = "o")
plot(b~x, xs.avg[[2]], type = "l", 
     main = "Reach average cross section 2",
     xlim = c(0,800), ylim = c(135,144))
plot(WSE~w, rWSEw[[2]], xlab = "width (m)", ylab = "height (m)", 
     ylim = c(135, 144), main = "Height-width relationship", type = "o")
plot(b~x, xs.avg[[3]], type = "l", 
     main = "Reach average cross section 3",
     xlim = c(0,800), ylim = c(135,144))
plot(WSE~w, rWSEw[[3]], xlab = "width (m)", ylab = "height (m)", 
     ylim = c(135, 144), main = "Height-width relationship", type = "o")

# ------------------------------------------------------------------------------------------------
# Main experiments - model fitting at different exposure levels, known measurements
# Should use errors in variables method, but this can be refined later.
# Will do this with uncertain measurements later, preferably without MC simulation

expo <- seq(0.05, 0.95, length.out = 19) # exposure levels
n_exp_levels <- length(expo)
nr <- length(rWSEw)

# Run the parallel computations (see computation_time_calculator.xls)

ncores <- detectCores()
# registerDoMC(cores = ncores - 1)
registerDoMC(cores = 3)

begin.time <- Sys.time()
WSEw_val <- foreach(r = 1:nr, .combine = c) %dopar% {observe_par(r)} # returns a dummy value; the point is to save .rds files
print(Sys.time() - begin.time)

begin.time <- Sys.time()
lfval <- foreach(r = 1:nr, .combine = c) %dopar% {fit_linear_par(r)}
print(Sys.time() - begin.time)

begin.time <- Sys.time()
sbval <- foreach(r = 1:nr, .combine = c) %dopar% {fit_sb_par(r)}
print(Sys.time() - begin.time)

begin.time <- Sys.time()
nlval <- foreach(r = 1:nr, .combine = c) %dopar% {fit_nl_par(r)}
print(Sys.time() - begin.time)

begin.time <- Sys.time() # sometimes the gradient is singular for the optimization routine in nlsLM
nlsbval <- foreach(r = 1:nr, .combine = c) %dopar% {fit_nlsb_par(r)}
print(Sys.time() - begin.time)

# --------------------------------------------------------------------------------------------------------------------------------------------
# Plot for presentation (one cross section)
r <- 3
WSEw_obs1 <- observe(rWSEw[[r]], exposure = 0.5, sd_w = 0, sd_wse = 0)
plot(WSE~w, WSEw_obs1, xlim = c(0, 800), ylim = c(131,143), 
     main = paste("Reach averaged cross section", r, "at 50% exposure"),
     xlab = "w (m)", ylab = "WSE (m)")
lines(WSE~w, rWSEw[[r]])
# points(0, rWSEw[[1]]$WSE[1], col = "black", pch = 8)

# Fit models
lf1 <- fit_linear(WSEw_obs1)
sb1 <- fit_slopebreak(WSEw_obs1, multiple_breaks = FALSE, continuity = TRUE)
nl1 <- fit_nonlinear(WSEw_obs1)
nlsb1 <- fit_nlsb(WSEw_obs1)

# plot linear
lines(WSEw_obs1$w, predict(lf1), col = "red")
w.vals <- seq(0, min(WSEw_obs1$w), 1)
lines(w.vals, predict(lf1, newdata = data.frame(w=w.vals)), 
      col = "red", lty = 2)
# points(0, predict(lf1, newdata = data.frame(w=0)), col = "red", pch = 19)

# plot slope break
nn <- length(WSEw_obs1$w)
sb1.ind <- attributes(sb1)$sb.ind
lines(WSEw_obs1$w[1:sb1.ind], predict(sb1[[1]]), col = "orange")
lines(WSEw_obs1$w[(sb1.ind):nn], predict(sb1[[2]]), col = "orange")
lines(w.vals, predict(sb1[[1]], newdata = data.frame(w=w.vals)), 
      col = "orange", lty = 2)
# points(0, predict(sb1[[1]], newdata = data.frame(w=0)), col = "orange", pch = 19)

# plot nl
lines(WSEw_obs1$w, predict(nl1), col = "green")
lines(w.vals, predict(nl1, newdata = data.frame(w=w.vals)), 
      col = "green", lty = 2)
# points(0, predict(nl1, newdata = data.frame(w=0)), col = "green", pch = 19)

# plot nlsb
nlsb1.ind <- attributes(nlsb1)$sb.ind
lines(WSEw_obs1$w[1:nlsb1.ind], predict(nlsb1[[1]]), col = "blue")
lines(WSEw_obs1$w[(nlsb1.ind):nn], predict(nlsb1[[2]]), col = "blue")
lines(w.vals, predict(nlsb1[[1]], newdata = data.frame(w=w.vals)), 
      col = "blue", lty = 2)
# points(0, predict(nlsb1[[1]], newdata = data.frame(w=0)), col = "blue", pch = 19)

legend("topleft", legend = c("Data", "Linear","SB","NL","NLSB"), 
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# ------------------------------------------------------------------------------------------------------------
# Compute true z0

z0.true.ra <- unlist(lapply(cross_sections_avg$b, min))
saveRDS(z0.true.ra, file = file.path(exp_dir, "z0_true_ra.rds")) # cross section

# ------------------------------------------------------------------------------------------------------------
# Compute true s0

# s0.true.xs <- diff(z0.true.xs)
# s0.true.ra <- diff(z0.true.ra)
# par(mfrow=c(2,1))
# plot(s0.true.xs, type = "l")
# plot(s0.true.ra, type = "l")
# abline(0,0)
# saveRDS(s0.true.xs, file = file.path(exp_dir, "s0_true_xs.rds")) 
# saveRDS(s0.true.ra, file = file.path(exp_dir, "s0_true_ra.rds"))

# ------------------------------------------------------------------------------------------------------------
# Compute true A, WP

# Calculate true hydraulic parameters for cross sections
# n.xs <- length(xWSEw)
# A.true.xs <- vector(length = n.xs)
# WP.true.xs <- vector(length = n.xs)
# for (k in 1:n.xs)
# {
#   x <- cross_sections$x[[k]]
#   b <- cross_sections$b[[k]]
#   WSE <- max(b) # bankfull
#   A.true.xs[k] <- calc_A(x, b, WSE)
#   WP.true.xs[k] <- calc_WP(x, b, WSE)
# }
# plot(A.true.xs, type="l")
# plot(WP.true.xs, type="l")
# saveRDS(A.true.xs, file = file.path(exp_dir, "A_true_xs.rds"))
# saveRDS(WP.true.xs, file = file.path(exp_dir, "WP_true_xs.rds"))

# Calculate reach-average hydraulic parameters
A.true.ra <- vector(length = nr)
WP.true.ra <- vector(length = nr)
for (k in 1:nr)
{
  x <- cross_sections_avg$x[[k]]
  b <- cross_sections_avg$b[[k]]
  WSE <- max(b) # bankfull
  A.true.ra[k] <- calc_A(x, b, WSE)
  WP.true.ra[k] <- calc_WP(x, b, WSE)
}
plot(A.true.ra, type="l")
plot(WP.true.ra, type="l")
saveRDS(A.true.ra, file = file.path(exp_dir, "A_true_ra.rds"))
saveRDS(WP.true.ra, file = file.path(exp_dir, "WP_true_ra.rds"))

# ------------------------------------------------------------------------------------------------------------
# Compute true A0

# Calculate true A0 (cross section)
# A0.true.xs <- array(dim = c(nr, n_exp_levels))
# w0.xs <- array(dim = c(nr, n_exp_levels)) # lowest observed width value
# for (r in 1:nr)
# {
#   for (k in 1:n_exp_levels)
#   {
#     WSEw_obs <- observe(xWSEw[[r]], sd_wse = 0, sd_w = 0, exposure = expo[k])
#     p <- max(which(xWSEw[[r]]$WSE<min(WSEw_obs$WSE))) # index up to which is not observed
#     w0.xs[r,k] <- min(WSEw_obs$w)
#     A0.true.xs[r,k] <- calc_A_from_WSEw(xWSEw[[r]][1:p,])
#   }
#   print(r)
# }
# saveRDS(w0.xs, file = file.path(exp_dir, "w0_xs.rds"))
# saveRDS(A0.true.xs, file = file.path(exp_dir, "A0_true_xs.rds"))

# Calculate true A0 (reach average)
A0.true.ra <- array(dim = c(nr, n_exp_levels))
w0.ra <- array(dim = c(nr, n_exp_levels)) # lowest observed width value
for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    WSEw_obs <- observe(rWSEw[[r]], sd_wse = 0, sd_w = 0, exposure = expo[k])
    p <- max(which(rWSEw[[r]]$WSE<min(WSEw_obs$WSE))) # index up to which is not observed
    w0.ra[r,k] <- min(WSEw_obs$w)
    A0.true.ra[r,k] <- calc_A_from_WSEw(rWSEw[[r]][1:p,])
  }
  print(r)
}
saveRDS(w0.ra, file = file.path(exp_dir, "w0_ra.rds"))
saveRDS(A0.true.ra, file = file.path(exp_dir, "A0_true_ra.rds"))

# ------------------------------------------------------------------------------------------------
# Predict hydraulic parameters

# Do computations in parallel
ncores <- detectCores()
registerDoMC(cores = ncores - 1)
registerDoMC(cores = 12)

begin.time <- Sys.time()
pred_lf <- foreach(r = 1:nr) %dopar% {pred_linear_par(r)}
print(Sys.time() - begin.time)
save(pred_lf, file = file.path(exp_dir, "pred_lf_bu.rda"))

begin.time <- Sys.time()
pred_sb <- foreach(r = 1:nr) %dopar% {pred_sb_par(r)}
print(Sys.time() - begin.time)
save(pred_sb, file = file.path(exp_dir, "pred_sb_bu.rda"))

begin.time <- Sys.time()
pred_nl <- foreach(r = 1:nr) %dopar% {pred_nl_par(r)}
print(Sys.time() - begin.time)
save(pred_nl, file = file.path(exp_dir, "pred_nl_bu.rda"))

begin.time <- Sys.time()
pred_nlsb <- foreach(r = 1:nr) %dopar% {pred_nlsb_par(r)}
print(Sys.time() - begin.time)
save(pred_nlsb, file = file.path(exp_dir, "pred_nlsb_bu.rda"))
# load(file.path(exp_dir, "pred_nlsb_bu.rda"))

# Run predictions again for any reaches where errors occurred

# ------------------------------------------------------------------------------------------------------------
# Reformat to a more convenient format

# Initialize
z0.l <- array(dim = c(nr, n_exp_levels, M))
z0.sb <- array(dim = c(nr, n_exp_levels, M))
z0.nl <- array(dim = c(nr, n_exp_levels, M))
z0.nlsb <- array(dim = c(nr, n_exp_levels, M))

A.l <- array(dim = c(nr, n_exp_levels, M))
A.sb <- array(dim = c(nr, n_exp_levels, M))
A.nl <- array(dim = c(nr, n_exp_levels, M))
A.nlsb <- array(dim = c(nr, n_exp_levels, M))

WP.l <- array(dim = c(nr, n_exp_levels, M))
WP.sb <- array(dim = c(nr, n_exp_levels, M))
WP.nl <- array(dim = c(nr, n_exp_levels, M))
WP.nlsb <- array(dim = c(nr, n_exp_levels, M))

A0.l <- array(dim = c(nr, n_exp_levels, M))
A0.sb <- array(dim = c(nr, n_exp_levels, M))
A0.nl <- array(dim = c(nr, n_exp_levels, M))
A0.nlsb <- array(dim = c(nr, n_exp_levels, M))

for (r in 1:nr)
{
  z0.l[r,,] <- pred_lf[[r]]$z0
  z0.sb[r,,] <- pred_sb[[r]]$z0
  z0.nl[r,,] <- pred_nl[[r]]$z0
  z0.nlsb[r,,] <- pred_nlsb[[r]]$z0
  
  A0.l[r,,] <- pred_lf[[r]]$A0
  A0.sb[r,,] <- pred_sb[[r]]$A0
  A0.nl[r,,] <- pred_nl[[r]]$A0
  A0.nlsb[r,,] <- pred_nlsb[[r]]$A0
  
  A.l[r,,] <- pred_lf[[r]]$A
  A.sb[r,,] <- pred_sb[[r]]$A
  A.nl[r,,] <- pred_nl[[r]]$A
  A.nlsb[r,,] <- pred_nlsb[[r]]$A
  
  WP.l[r,,] <- pred_lf[[r]]$WP
  WP.sb[r,,] <- pred_sb[[r]]$WP
  WP.nl[r,,] <- pred_nl[[r]]$WP
  WP.nlsb[r,,] <- pred_nlsb[[r]]$WP
}

save(z0.l, z0.sb, z0.nl, z0.nlsb, file = file.path(exp_dir, "z0_pred.rda"))
save(A.l, A.sb, A.nl, A.nlsb, file = file.path(exp_dir, "A_pred.rda"))
save(WP.l, WP.sb, WP.nl, WP.nlsb, file = file.path(exp_dir, "WP_pred.rda"))
save(A0.l, A0.sb, A0.nl, A0.nlsb, file = file.path(exp_dir, "A0_pred.rda"))

# Compute slope via finite difference
# s0.l <- apply(z0.l, c(2,3), diff)
# s0.sb <- apply(z0.sb, c(2,3), diff)
# s0.nl <- apply(z0.nl, c(2,3), diff)
# s0.nlsb <- apply(z0.nlsb, c(2,3), diff)
# save(s0.l, s0.sb, s0.nl, s0.nlsb, file = file.path(exp_dir, "s0_pred.rda"))

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ---------------------------------------- Plotting ----------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Consider switching to Matlab and using tools from CEE251D DA project - nah, better to stick to one language

par(mfrow = c(1,1), mar = c(5,5,5,5))
par(mfrow = c(2,2))

# ------------------------------------------------------------------------------------------------
# Test for Gaussianity

# Batch test
r <- 1
k <- 19
z <- apply(s0.nlsb, c(1,2), is_normal)
z
hist(s0.nlsb[r, k,], "fd")

# The linear model z0 predictions appear Gaussian, but not the others, use quantiles in plots.

# ------------------------------------------------------------------------------------------------
# Estimate moments of sampling distributions of predicted values

# z0
z0.l.lower <- array(dim = c(nr, n_exp_levels))
z0.l.med <- array(dim = c(nr, n_exp_levels))
z0.l.upper <- array(dim = c(nr, n_exp_levels))

z0.sb.lower <- array(dim = c(nr, n_exp_levels))
z0.sb.med <- array(dim = c(nr, n_exp_levels))
z0.sb.upper <- array(dim = c(nr, n_exp_levels))

z0.nl.lower <- array(dim = c(nr, n_exp_levels))
z0.nl.med <- array(dim = c(nr, n_exp_levels))
z0.nl.upper <- array(dim = c(nr, n_exp_levels))

z0.nlsb.lower <- array(dim = c(nr, n_exp_levels))
z0.nlsb.med <- array(dim = c(nr, n_exp_levels))
z0.nlsb.upper <- array(dim = c(nr, n_exp_levels))

for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    
    q1 <- quantile(z0.l[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    z0.l.lower[r,k] <- as.numeric(q1[1])
    z0.l.med[r,k] <- as.numeric(q1[2])
    z0.l.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(z0.sb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    z0.sb.lower[r,k] <- as.numeric(q1[1])
    z0.sb.med[r,k] <- as.numeric(q1[2])
    z0.sb.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(z0.nl[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    z0.nl.lower[r,k] <- as.numeric(q1[1])
    z0.nl.med[r,k] <- as.numeric(q1[2])
    z0.nl.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(z0.nlsb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    z0.nlsb.lower[r,k] <- as.numeric(q1[1])
    z0.nlsb.med[r,k] <- as.numeric(q1[2])
    z0.nlsb.upper[r,k] <- as.numeric(q1[3])
  }
}

# A0
A0.l.lower <- array(dim = c(nr, n_exp_levels))
A0.l.med <- array(dim = c(nr, n_exp_levels))
A0.l.upper <- array(dim = c(nr, n_exp_levels))

A0.sb.lower <- array(dim = c(nr, n_exp_levels))
A0.sb.med <- array(dim = c(nr, n_exp_levels))
A0.sb.upper <- array(dim = c(nr, n_exp_levels))

A0.nl.lower <- array(dim = c(nr, n_exp_levels))
A0.nl.med <- array(dim = c(nr, n_exp_levels))
A0.nl.upper <- array(dim = c(nr, n_exp_levels))

A0.nlsb.lower <- array(dim = c(nr, n_exp_levels))
A0.nlsb.med <- array(dim = c(nr, n_exp_levels))
A0.nlsb.upper <- array(dim = c(nr, n_exp_levels))

for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    q1 <- quantile(A0.l[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A0.l.lower[r,k] <- as.numeric(q1[1])
    A0.l.med[r,k] <- as.numeric(q1[2])
    A0.l.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A0.sb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A0.sb.lower[r,k] <- as.numeric(q1[1])
    A0.sb.med[r,k] <- as.numeric(q1[2])
    A0.sb.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A0.nl[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A0.nl.lower[r,k] <- as.numeric(q1[1])
    A0.nl.med[r,k] <- as.numeric(q1[2])
    A0.nl.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A0.nlsb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A0.nlsb.lower[r,k] <- as.numeric(q1[1])
    A0.nlsb.med[r,k] <- as.numeric(q1[2])
    A0.nlsb.upper[r,k] <- as.numeric(q1[3])
  }
}

# s0
s0.l.lower <- array(dim = c(nr, n_exp_levels))
s0.l.med <- array(dim = c(nr, n_exp_levels))
s0.l.upper <- array(dim = c(nr, n_exp_levels))

s0.sb.lower <- array(dim = c(nr, n_exp_levels))
s0.sb.med <- array(dim = c(nr, n_exp_levels))
s0.sb.upper <- array(dim = c(nr, n_exp_levels))

s0.nl.lower <- array(dim = c(nr, n_exp_levels))
s0.nl.med <- array(dim = c(nr, n_exp_levels))
s0.nl.upper <- array(dim = c(nr, n_exp_levels))

s0.nlsb.lower <- array(dim = c(nr, n_exp_levels))
s0.nlsb.med <- array(dim = c(nr, n_exp_levels))
s0.nlsb.upper <- array(dim = c(nr, n_exp_levels))

for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    q1 <- quantile(s0.l[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    s0.l.lower[r,k] <- as.numeric(q1[1])
    s0.l.med[r,k] <- as.numeric(q1[2])
    s0.l.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(s0.sb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    s0.sb.lower[r,k] <- as.numeric(q1[1])
    s0.sb.med[r,k] <- as.numeric(q1[2])
    s0.sb.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(s0.nl[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    s0.nl.lower[r,k] <- as.numeric(q1[1])
    s0.nl.med[r,k] <- as.numeric(q1[2])
    s0.nl.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(s0.nlsb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    s0.nlsb.lower[r,k] <- as.numeric(q1[1])
    s0.nlsb.med[r,k] <- as.numeric(q1[2])
    s0.nlsb.upper[r,k] <- as.numeric(q1[3])
  }
}

# WP
WP.l.lower <- array(dim = c(nr, n_exp_levels))
WP.l.med <- array(dim = c(nr, n_exp_levels))
WP.l.upper <- array(dim = c(nr, n_exp_levels))

WP.sb.lower <- array(dim = c(nr, n_exp_levels))
WP.sb.med <- array(dim = c(nr, n_exp_levels))
WP.sb.upper <- array(dim = c(nr, n_exp_levels))

WP.nl.lower <- array(dim = c(nr, n_exp_levels))
WP.nl.med <- array(dim = c(nr, n_exp_levels))
WP.nl.upper <- array(dim = c(nr, n_exp_levels))

WP.nlsb.lower <- array(dim = c(nr, n_exp_levels))
WP.nlsb.med <- array(dim = c(nr, n_exp_levels))
WP.nlsb.upper <- array(dim = c(nr, n_exp_levels))

for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    q1 <- quantile(WP.l[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    WP.l.lower[r,k] <- as.numeric(q1[1])
    WP.l.med[r,k] <- as.numeric(q1[2])
    WP.l.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(WP.sb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    WP.sb.lower[r,k] <- as.numeric(q1[1])
    WP.sb.med[r,k] <- as.numeric(q1[2])
    WP.sb.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(WP.nl[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    WP.nl.lower[r,k] <- as.numeric(q1[1])
    WP.nl.med[r,k] <- as.numeric(q1[2])
    WP.nl.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(WP.nlsb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    WP.nlsb.lower[r,k] <- as.numeric(q1[1])
    WP.nlsb.med[r,k] <- as.numeric(q1[2])
    WP.nlsb.upper[r,k] <- as.numeric(q1[3])
  }
}

# A
A.l.lower <- array(dim = c(nr, n_exp_levels))
A.l.med <- array(dim = c(nr, n_exp_levels))
A.l.upper <- array(dim = c(nr, n_exp_levels))

A.sb.lower <- array(dim = c(nr, n_exp_levels))
A.sb.med <- array(dim = c(nr, n_exp_levels))
A.sb.upper <- array(dim = c(nr, n_exp_levels))

A.nl.lower <- array(dim = c(nr, n_exp_levels))
A.nl.med <- array(dim = c(nr, n_exp_levels))
A.nl.upper <- array(dim = c(nr, n_exp_levels))

A.nlsb.lower <- array(dim = c(nr, n_exp_levels))
A.nlsb.med <- array(dim = c(nr, n_exp_levels))
A.nlsb.upper <- array(dim = c(nr, n_exp_levels))

for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    q1 <- quantile(A.l[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A.l.lower[r,k] <- as.numeric(q1[1])
    A.l.med[r,k] <- as.numeric(q1[2])
    A.l.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A.sb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A.sb.lower[r,k] <- as.numeric(q1[1])
    A.sb.med[r,k] <- as.numeric(q1[2])
    A.sb.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A.nl[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A.nl.lower[r,k] <- as.numeric(q1[1])
    A.nl.med[r,k] <- as.numeric(q1[2])
    A.nl.upper[r,k] <- as.numeric(q1[3])
    
    q1 <- quantile(A.nlsb[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A.nlsb.lower[r,k] <- as.numeric(q1[1])
    A.nlsb.med[r,k] <- as.numeric(q1[2])
    A.nlsb.upper[r,k] <- as.numeric(q1[3])
  }
}

# |||||||||||||
# vvvvvvvvvvvvv
# Entered commands through here, starting 9/24 at 8:52 p.m.

# ------------------------------------------------------------------------------------------------
# Make plots of median z0, A0 along the river

par(mfrow=c(2,2))
k <- 16 # exposure level
l.pred <- data.frame(r = 1:nr, 
                     z0.med = z0.l.med[,k], z0.lower = z0.l.lower[,k], z0.upper = z0.l.upper[,k],
                     A0.med = A0.l.med[,k], A0.lower = A0.l.lower[,k], A0.upper = A0.l.upper[,k],
                     s0.med = s0.l.med[,k], s0.lower = s0.l.lower[,k], s0.upper = s0.l.upper[,k],
                     WP.med = WP.l.med[,k], WP.lower = WP.l.lower[,k], WP.upper = WP.l.upper[,k],
                     A.med = A.l.med[,k], A.lower = A.l.lower[,k], A.upper = A.l.upper[,k])
sb.pred <- data.frame(r = 1:nr, 
                      z0.med = z0.sb.med[,k], z0.lower = z0.sb.lower[,k], z0.upper = z0.sb.upper[,k],
                      A0.med = A0.sb.med[,k], A0.lower = A0.sb.lower[,k], A0.upper = A0.sb.upper[,k],
                      s0.med = s0.sb.med[,k], s0.lower = s0.sb.lower[,k], s0.upper = s0.sb.upper[,k],
                      WP.med = WP.sb.med[,k], WP.lower = WP.sb.lower[,k], WP.upper = WP.sb.upper[,k],
                      A.med = A.sb.med[,k], A.lower = A.sb.lower[,k], A.upper = A.sb.upper[,k])
nl.pred <- data.frame(r = 1:nr, 
                      z0.med = z0.nl.med[,k], z0.lower = z0.nl.lower[,k], z0.upper = z0.nl.upper[,k],
                      A0.med = A0.nl.med[,k], A0.lower = A0.nl.lower[,k], A0.upper = A0.nl.upper[,k],
                      s0.med = s0.nl.med[,k], s0.lower = s0.nl.lower[,k], s0.upper = s0.nl.upper[,k],
                      WP.med = WP.nl.med[,k], WP.lower = WP.nl.lower[,k], WP.upper = WP.nl.upper[,k],
                      A.med = A.nl.med[,k], A.lower = A.nl.lower[,k], A.upper = A.nl.upper[,k])
nlsb.pred <- data.frame(r = 1:nr, 
                        z0.med = z0.nlsb.med[,k], z0.lower = z0.nlsb.lower[,k], z0.upper = z0.nlsb.upper[,k],
                        A0.med = A0.nlsb.med[,k], A0.lower = A0.nlsb.lower[,k], A0.upper = A0.nlsb.upper[,k],
                        s0.med = s0.nlsb.med[,k], s0.lower = s0.nlsb.lower[,k], s0.upper = s0.nlsb.upper[,k],
                        WP.med = WP.nlsb.med[,k], WP.lower = WP.nlsb.lower[,k], WP.upper = WP.nlsb.upper[,k],
                        A.med = A.nlsb.med[,k], A.lower = A.nlsb.lower[,k], A.upper = A.nlsb.upper[,k])

# Do reach-averaging
l.pred.ra <- lapply(l.pred, ra, n = 2000)
sb.pred.ra <- lapply(sb.pred, ra, n = 2000)
nl.pred.ra <- lapply(nl.pred, ra, n = 2000)
nlsb.pred.ra <- lapply(nlsb.pred, ra, n = 2000)
save(l.pred.ra, sb.pred.ra, nl.pred.ra, nlsb.pred.ra, file = file.path(exp_dir, "pred_ra_80.rda")) # save predictions (ra)

# Somewhere, this is messing up the nonlinear predictions


load(file.path(exp_dir, "pred_ra_20.rda"))
load(file.path(exp_dir, "pred_ra_40.rda"))
load(file.path(exp_dir, "pred_ra_60.rda"))
load(file.path(exp_dir, "pred_ra_80.rda"))

# z0
plot(z0.true.ra, main = paste("Median z0 (m) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(100,150), lwd = 1.5,
     xlab = "reach-average cross section", ylab = "z0")
lines(l.pred.ra$z0.med, col = "red")
lines(sb.pred.ra$z0.med, col = "orange")
lines(nl.pred.ra$z0.med, col = "green")
lines(nlsb.pred.ra$z0.med, col = "blue")
legend("bottomright", legend = c("True", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A0
plot(A0.true.ra[,k], main = paste("Median A0 at", expo[k]*100,"% exposure (sq. m)"), 
     type = "l", ylim = c(0,4500), lwd = 1.5,
     xlab = "Cross section", ylab = "A0")
lines(l.pred.ra$A0.med, col = "red")
lines(sb.pred.ra$A0.med, col = "orange")
lines(nl.pred.ra$A0.med, col = "green")
lines(nlsb.pred.ra$A0.med, col = "blue")
legend("bottomright", legend = c("True", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# s0
plot(s0.true.ra*1e5, main = paste("s0 (cm/km) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(-2000,2000), 
     xlab = "reach-average cross section", ylab = "bed slope")
lines(1e5*l.pred$s0.med, col = "red")
lines(1e5*sb.pred$s0.med, col = "orange")
lines(1e5*nl.pred$s0.med, col = "green")
lines(1e5*nlsb.pred$s0.med, col = "blue")
legend("topleft", legend = c("True", "Linear","SB","NL","NLSB"), 
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

K <- 1000 # Smooth slope using a k-point moving average
s0.true.ra.smooth <- filter(s0.true.ra, sides = 2, filter = rep(1/K,K))
s0.l.smooth <- filter(l.pred$s0.med, sides = 2, filter = rep(1/K,K))
s0.sb.smooth <- filter(sb.pred$s0.med, sides = 2, filter = rep(1/K,K))
s0.nl.smooth <- filter(nl.pred$s0.med, sides = 2, filter = rep(1/K,K))
s0.nlsb.smooth <- filter(nlsb.pred$s0.med, sides = 2, filter = rep(1/K,K))

# s0 smoothed
plot(s0.true.ra.smooth*1e5, main = paste("s0 (cm/km) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(-500,500), 
     xlab = "reach-average cross section", ylab = "smoothed bed slope")
lines(1e5*s0.l.smooth, col = "red")
lines(1e5*s0.sb.smooth, col = "orange")
lines(1e5*s0.nl.smooth, col = "green")
lines(1e5*s0.nlsb.smooth, col = "blue")
legend("bottomright", legend = c("True", "Linear","SB","NL","NLSB"), 
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A
plot(A.true.ra[1:10], main = paste("Flow area (m^2) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(0,6000))
lines(l.pred$A.med, col = "red")
lines(sb.pred$A.med, col = "orange")
lines(nl.pred$A.med, col = "green")
lines(nlsb.pred$A.med, col = "blue")
legend("bottomleft", legend = c("True", "Linear","SB","NL","NLSB"), 
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# WP
plot(WP.true.ra[1:10], main = paste("Wetted perimeter (m) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(450,850))
lines(l.pred$WP.med, col = "red")
lines(sb.pred$WP.med, col = "orange")
lines(nl.pred$WP.med, col = "green")
lines(nlsb.pred$WP.med, col = "blue")
legend("topleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# ------------------------------------------------------------------------------------------------
# Make plots of average z0, A, WP, s0, A0 error at each exposure level

# z0
z0.bias <- plot_bias(expo, z0_med, z0.true.ra[1:nr], na.rm = TRUE,
                     main = "z0 bias vs. exposure level, no meas. error", ylab = "Bias (m)", ylim = c(-25,3))
legend("bottomright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# s0
s0.bias <- plot_bias(expo, s0_med, s0.true.ra[1:(nr-1)], na.rm = TRUE,
                     main = "s0 bias vs. exposure level, no meas. error", ylab = "Bias (cm/km)")
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A
A.bias <- plot_bias(expo, A_med, A.true.ra[1:nr], na.rm = TRUE,
                    main = "A bias vs. exposure level, no meas. error", ylab = "Bias (m^2)")
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# WP
WP.bias <- plot_bias(expo, WP_med, WP.true.ra[1:nr], na.rm = TRUE,
                     main = "WP bias vs. exposure level, no meas. error", ylab = "Bias (m)")
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A0
A0.bias <- array(dim = c(n_exp_levels, 4)) # cannot use plot_bias because A0 changes with exposure level
na.rm = TRUE
for (k in 1:n_exp_levels)
{
  A0.bias[k,1] <- mean(A0.l.med - A0.true.ra[,k], na.rm = na.rm)
  A0.bias[k,2] <- mean(A0.sb.med - A0.true.ra[,k], na.rm = na.rm)
  A0.bias[k,3] <- mean(A0.nl.med - A0.true.ra[,k], na.rm = na.rm)
  A0.bias[k,4] <- mean(A0.nlsb.med - A0.true.ra[,k], na.rm = na.rm)
}
A0.bias <- as.data.frame(A0.bias)
names(A0.bias) <- c("l","sb","nl","nlsb")

# Plot A0.bias vs. exposure level
par(opar)
plot(100*expo, abs(A0.bias$l), col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "Abs(A0.bias) using median as estimate (m2)", ylab = "A0.bias (m)", ylim = c(0,2000))
lines(100*expo, abs(A0.bias$sb), col = "orange")
lines(100*expo, abs(A0.bias$nl), col = "green")
lines(100*expo, abs(A0.bias$nlsb), col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

save(z0.bias, s0.bias, A.bias, WP.bias, A0.bias, file = file.path(exp_dir, "bias.rda"))

# ------------------------------------------------------------------------------------------------
# Compute RMSE

z0.rmse <- compute_rmse(z0_med, z0.true.ra)
A.rmse <- compute_rmse(A_med, A.true.ra)
WP.rmse <- compute_rmse(WP_med, WP.true.ra)
s0.rmse <- compute_rmse(s0_med, s0.true.ra)

A0.rmse <- array(dim = c(n_exp_levels, 4)) # cannot use compute_rmse because A0 changes with exposure level
na.rm = TRUE
for (k in 1:n_exp_levels)
{
  A0.rmse[k,1] <- ((1/length(A0.true.ra[,k]))*t(A0.l.med[,k] - A0.true.ra[,k])%*%(A0.l.med[,k] - A0.true.ra[,k]))^0.5
  A0.rmse[k,2] <- ((1/length(A0.true.ra[,k]))*t(A0.sb.med[,k] - A0.true.ra[,k])%*%(A0.sb.med[,k] - A0.true.ra[,k]))^0.5
  A0.rmse[k,3] <- ((1/length(A0.true.ra[,k]))*t(A0.nl.med[,k] - A0.true.ra[,k])%*%(A0.nl.med[,k] - A0.true.ra[,k]))^0.5
  A0.rmse[k,4] <- ((1/length(A0.true.ra[,k]))*t(A0.nlsb.med[,k] - A0.true.ra[,k])%*%(A0.nlsb.med[,k] - A0.true.ra[,k]))^0.5
}
A0.rmse <- as.data.frame(A0.rmse)
names(A0.rmse) <- c("l","sb","nl","nlsb")

save(z0.rmse, s0.rmse, A.rmse, WP.rmse, A0.rmse, file = file.path(exp_dir, "rmse.rda"))

# ------------------------------------------------------------------------------------------------
# Make plots of z0, A, WP, s0, A0 RMSE at each exposure level

# z0
plot(100*expo, z0.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "z0 RMSE (m)", ylab = "z0.rmse (m)",
     ylim = c(0,10))
lines(100*expo, z0.rmse$sb, col = "orange")
lines(100*expo, z0.rmse$nl, col = "green")
lines(100*expo, z0.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# s0
plot(100*expo, s0.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "s0 RMSE (cm/km)", ylab = "s0.rmse (m)")
lines(100*expo, s0.rmse$sb, col = "orange")
lines(100*expo, s0.rmse$nl, col = "green")
lines(100*expo, s0.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A
plot(100*expo, A.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "A RMSE (m2)", ylab = "A.rmse (m)",
     ylim = c(0,7000))
lines(100*expo, A.rmse$sb, col = "orange")
lines(100*expo, A.rmse$nl, col = "green")
lines(100*expo, A.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# WP
plot(100*expo, WP.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "WP RMSE (m)", ylab = "WP.rmse (m)",
     ylim = c(0,4))
lines(100*expo, WP.rmse$sb, col = "orange")
lines(100*expo, WP.rmse$nl, col = "green")
lines(100*expo, WP.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# A0
plot(100*expo, A0.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)",
     main = "A0.rmse (m2)", ylab = "A0.rmse (m)", ylim = c(0,2000))
lines(100*expo, A0.rmse$sb, col = "orange")
lines(100*expo, A0.rmse$nl, col = "green")
lines(100*expo, A0.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","NL","NLSB"),
       col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)

# ------------------------------------------------------------------------------------------------
# Make plots of median z0, A0 profiles with error bounds

par(mfrow = c(2,2))

k <- 16 # exposure level
l.pred <- data.frame(r = 1:nr, 
                     z0.med = z0.l.med[,k], z0.lower = z0.l.lower[,k], z0.upper = z0.l.upper[,k],
                     A0.med = A0.l.med[,k], A0.lower = A0.l.lower[,k], A0.upper = A0.l.upper[,k],
                     s0.med = s0.l.med[,k], s0.lower = s0.l.lower[,k], s0.upper = s0.l.upper[,k],
                     WP.med = WP.l.med[,k], WP.lower = WP.l.lower[,k], WP.upper = WP.l.upper[,k],
                     A.med = A.l.med[,k], A.lower = A.l.lower[,k], A.upper = A.l.upper[,k])
sb.pred <- data.frame(r = 1:nr, 
                      z0.med = z0.sb.med[,k], z0.lower = z0.sb.lower[,k], z0.upper = z0.sb.upper[,k],
                      A0.med = A0.sb.med[,k], A0.lower = A0.sb.lower[,k], A0.upper = A0.sb.upper[,k],
                      s0.med = s0.sb.med[,k], s0.lower = s0.sb.lower[,k], s0.upper = s0.sb.upper[,k],
                      WP.med = WP.sb.med[,k], WP.lower = WP.sb.lower[,k], WP.upper = WP.sb.upper[,k],
                      A.med = A.sb.med[,k], A.lower = A.sb.lower[,k], A.upper = A.sb.upper[,k])
nl.pred <- data.frame(r = 1:nr, 
                      z0.med = z0.nl.med[,k], z0.lower = z0.nl.lower[,k], z0.upper = z0.nl.upper[,k],
                      A0.med = A0.nl.med[,k], A0.lower = A0.nl.lower[,k], A0.upper = A0.nl.upper[,k],
                      s0.med = s0.nl.med[,k], s0.lower = s0.nl.lower[,k], s0.upper = s0.nl.upper[,k],
                      WP.med = WP.nl.med[,k], WP.lower = WP.nl.lower[,k], WP.upper = WP.nl.upper[,k],
                      A.med = A.nl.med[,k], A.lower = A.nl.lower[,k], A.upper = A.nl.upper[,k])
nlsb.pred <- data.frame(r = 1:nr, 
                        z0.med = z0.nlsb.med[,k], z0.lower = z0.nlsb.lower[,k], z0.upper = z0.nlsb.upper[,k],
                        A0.med = A0.nlsb.med[,k], A0.lower = A0.nlsb.lower[,k], A0.upper = A0.nlsb.upper[,k],
                        s0.med = s0.nlsb.med[,k], s0.lower = s0.nlsb.lower[,k], s0.upper = s0.nlsb.upper[,k],
                        WP.med = WP.nlsb.med[,k], WP.lower = WP.nlsb.lower[,k], WP.upper = WP.nlsb.upper[,k],
                        A.med = A.nlsb.med[,k], A.lower = A.nlsb.lower[,k], A.upper = A.nlsb.upper[,k])

# z0
plot(z0.true.ra, main = paste("Median z0 (m) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(120,140), lwd = 1.5,
     xlab = "reach-average cross section", ylab = "z0")
lines(nl.pred$z0.lower, col = "lightgreen", lty = 2)
lines(nl.pred$z0.med, col = "green")
lines(nl.pred$z0.upper, col = "lightgreen", lty = 2)
lines(nlsb.pred$z0.lower, col = "lightblue", lty = 2)
lines(nlsb.pred$z0.med, col = "blue")
lines(nlsb.pred$z0.upper, col = "lightblue", lty = 2)
legend("bottomright", legend = c("True","NL","NLSB"),
       col = c("black","green","blue"), lwd = c(1,1,1))

# A0
plot(A0.true.ra[,k], main = paste("Median A0 at", expo[k]*100,"% exposure (sq. m)"), 
     type = "l", ylim = c(0,500), lwd = 1.5,
     xlab = "Cross section", ylab = "A0")
lines(nl.pred$A0.lower, col = "lightgreen", lty = 2)
lines(nl.pred$A0.med, col = "green")
lines(nl.pred$A0.upper, col = "lightgreen", lty = 2)
lines(nlsb.pred$A0.lower, col = "lightblue", lty = 2)
lines(nlsb.pred$A0.med, col = "blue")
lines(nlsb.pred$A0.upper, col = "lightblue", lty = 2)
legend("bottomright", legend = c("True","NL","NLSB"),
       col = c("black", "green","blue"), lwd = c(1,1,1))

# -----------------------------------------------------------------------------------------------------------------------------------
# Calculate average spread for A0 and z0 predictions

k <- 16
# mean(z0.l.upper[,k] - z0.l.lower[,k], na.rm = TRUE)
# mean(z0.sb.upper[,k] - z0.sb.lower[,k], na.rm = TRUE)
# # mean(z0.sbm.upper[,k] - z0.sbm.lower[,k], na.rm = TRUE)
# mean(z0.nl.upper[,k] - z0.nl.lower[,k], na.rm = TRUE)
# mean(z0.nlsb.upper[,k] - z0.nlsb.lower[,k], na.rm = TRUE)

# z0 IQR
z0.l.iqr <- apply(z0.l, c(1,2), IQR, na.rm = TRUE)
z0.sb.iqr <- apply(z0.sb, c(1,2), IQR, na.rm = TRUE)
z0.nl.iqr <- apply(z0.nl, c(1,2), IQR, na.rm = TRUE)
z0.nlsb.iqr <- apply(z0.nlsb, c(1,2), IQR, na.rm = TRUE)

# A0 IQR
A0.l.iqr <- apply(A0.l, c(1,2), IQR, na.rm = TRUE)
A0.sb.iqr <- apply(A0.sb, c(1,2), IQR, na.rm = TRUE)
A0.nl.iqr <- apply(A0.nl, c(1,2), IQR, na.rm = TRUE)
A0.nlsb.iqr <- apply(A0.nlsb, c(1,2), IQR, na.rm = TRUE)

# z0 sd
z0.l.sd <- apply(z0.l, c(1,2), sd, na.rm = TRUE)
z0.sb.sd <- apply(z0.sb, c(1,2), sd, na.rm = TRUE)
z0.nl.sd <- apply(z0.nl, c(1,2), sd, na.rm = TRUE)
z0.nlsb.sd <- apply(z0.nlsb, c(1,2), sd, na.rm = TRUE)

# A0 sd
A0.l.sd <- apply(A0.l, c(1,2), sd, na.rm = TRUE)
A0.sb.sd <- apply(A0.sb, c(1,2), sd, na.rm = TRUE)
A0.nl.sd <- apply(A0.nl, c(1,2), sd, na.rm = TRUE)
A0.nlsb.sd <- apply(A0.nlsb, c(1,2), sd, na.rm = TRUE)

# calculate mean IQR and sd for several different exposure levels
meanIQR <- vector(length = 4, "list") # first do for z0, then for A0
meanSD <- vector(length = 4, "list")
names(meanIQR) <- c("l","sb","nl","nlsb")
names(meanSD) <- c("l","sb","nl","nlsb")
for (k in 1:n_exp_levels)
{
  meanIQR$l[k] <- mean(A0.l.iqr[,k], na.rm = TRUE)
  meanIQR$sb[k] <- mean(A0.sb.iqr[,k], na.rm = TRUE)
  meanIQR$nl[k] <- mean(A0.nl.iqr[,k], na.rm = TRUE)
  meanIQR$nlsb[k] <- mean(A0.nlsb.iqr[,k], na.rm = TRUE)
  
  meanSD$l[k] <- mean(A0.l.sd[,k], na.rm = TRUE)
  meanSD$sb[k] <- mean(A0.sb.sd[,k], na.rm = TRUE)
  meanSD$nl[k] <- mean(A0.nl.sd[,k], na.rm = TRUE)
  meanSD$nlsb[k] <- mean(A0.nlsb.sd[,k], na.rm = TRUE)
}

# Check for normality: the IQR should be 1.34898*SD for a normally distributed RV
326*1.34898 # many of the predicted hydraulic variables appear Gaussian, based on this test.
# Should also do visual inspection

hist()
hist(A0.l[,,8])
hist(A0.nlsb[,16,])
IQR(A0.nlsb[,16,], na.rm = TRUE)
sd(A0.nlsb[,16,], na.rm = TRUE)

# -----------------------------------------------------------------------------------------------------------------------------------
# Plot prediction error histograms


# tmp <- 
# hist(tmp[,k,], main = paste("z0 error (m) at", expo[k]*100,"% exposure"), xlab = "z0 prediction error")
# sd(tmp[,k,])
# IQR(tmp[,k,])

# z0 prediction error
par(mfrow = c(4,4))
k <- 16
xl <- -4
xu <- 2
hist((z0.l - z0.true.ra)[,k,], 
     main = paste("z0 error (m) at", expo[k]*100,"% exposure, linear method"), 
     xlab = "z0 prediction error", "fd", xlim = c(xl,xu))
hist((z0.sb - z0.true.ra)[,k,], 
     main = paste("z0 error (m) at", expo[k]*100,"% exposure, SB method"), 
     xlab = "z0 prediction error", "fd", xlim = c(xl,xu))
hist((z0.nl - z0.true.ra)[,k,], 
     main = paste("z0 error (m) at", expo[k]*100,"% exposure, NL method"), 
     xlab = "z0 prediction error", "fd", xlim = c(xl,xu))
hist((z0.nlsb - z0.true.ra)[,k,], 
     main = paste("z0 error (m) at", expo[k]*100,"% exposure, NLSB method"), 
     xlab = "z0 prediction error", "fd", xlim = c(xl,xu))

# A0 prediction error
par(mfrow = c(4,4))
k <- 16
xl <- -100
xu <- 400
hist((A0.l - A0.true.ra[,k])[,k,], 
     main = paste("A0 error (sq. m) at", expo[k]*100,"% exposure, linear method"), 
     xlab = "A0 prediction error", "fd", xlim = c(xl,xu))
hist((A0.sb - A0.true.ra[,k])[,k,], 
     main = paste("A0 error (sq. m) at", expo[k]*100,"% exposure, SB method"), 
     xlab = "A0 prediction error", "fd", xlim = c(xl,xu))
hist((A0.nl - A0.true.ra[,k])[,k,], 
     main = paste("A0 error (sq. m) at", expo[k]*100,"% exposure, NL method"), 
     xlab = "A0 prediction error", "fd", xlim = c(xl,xu))
hist((A0.nlsb - A0.true.ra[,k])[,k,], 
     main = paste("A0 error (sq. m) at", expo[k]*100,"% exposure, NLSB method"), 
     xlab = "A0 prediction error", "fd", xlim = c(xl,xu))

# ------------------------------------------------------------------------------------------
# Plot predicted z0 and A0 from each model for each reach-average cross section

par(mfrow = c(2,2))

for (k in seq(4,16,by=4))
{
  plot(xlim = c(0,4), 
       ylim = c(100, 144), xlab = "cross section", ylab = "z0 (m)",
       main = paste("z0 predictions,", 100*expo[k], "percent exposure"), 
       "n")
  if (k==16)
  {
    legend("bottomright", legend = c("Truth", "Linear","SB","NL","NLSB"),
           col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)
  }
  for (r in 1:nr)
  {
    rx <- seq(r-0.25, r+0.25, by = 0.01)
    ry <- rep(z0.true.ra[r], length(rx))
    lines(rx, ry, col = "black", pch = 7)
    points(r, z0.l[r,k,1], col = "red", pch = 1, cex = 1.5)
    points(r, z0.sb[r,k,1], col = "orange", pch = 1, cex = 1.5)
    points(r, z0.nl[r,k,1], col = "green", pch = 1, cex = 1.5)
    points(r, z0.nlsb[r,k,1], col = "blue", pch = 1, cex = 1.5)
  }
}

for (k in seq(4,16,by=4)) # something is wrong with the code here
{
  plot(xlim = c(0,4), 
       ylim = c(0, 4500), xlab = "cross section", ylab = "A0 (sq. m)",
       main = paste("A0 predictions,", 100*expo[k], "percent exposure"), 
       "n")
  if (k==16)
  {
    legend("topright", legend = c("Truth", "Linear","SB","NL","NLSB"),
           col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)
  }
  for (r in 1:nr)
  {
    rx <- seq(r-0.25, r+0.25, by = 0.01)
    ry <- rep(A0.true.ra[r], length(rx))
    lines(rx, ry, col = "black", pch = 7)
    points(r, A0.l[r,k,1], col = "red", pch = 1, cex = 1.5)
    points(r, A0.sb[r,k,1], col = "orange", pch = 1, cex = 1.5)
    points(r, A0.nl[r,k,1], col = "green", pch = 1, cex = 1.5)
    points(r, A0.nlsb[r,k,1], col = "blue", pch = 1, cex = 1.5)
  }
}

# ------------------------------------------------------------------------------------------------
# SCRAP 

# it would be good to take off outlier cross sections - say, those with bankfull parameters that are outliers
ra.xs <- vector(length = nr, "list")
# make this efficient later
sumw <- 0
sumWSE <- 0 
for (r in 1:2000)
{
  sumw <- sumw + xWSEw[[r]]$w
  sumWSE <- sumWSE + xWSEw[[r]]$WSE
}
ra.xs[[1]] <- data.frame(WSE = sumWSE/2000, w = sumw/2000)

sumw <- 0
sumWSE <- 0 
for (r in 2001:4000)
{
  sumw <- sumw + xWSEw[[r]]$w
  sumWSE <- sumWSE + xWSEw[[r]]$WSE
}
ra.xs[[2]] <- data.frame(WSE = sumWSE/2000, w = sumw/2000)

sumw <- 0
sumWSE <- 0 
for (r in 4001:5773)
{
  sumw <- sumw + xWSEw[[r]]$w
  sumWSE <- sumWSE + xWSEw[[r]]$WSE
}
# 5773-4001+1
ra.xs[[3]] <- data.frame(WSE = sumWSE/1773, w = sumw/1773)

plot(WSE~w, ra.xs[[1]])
plot(WSE~w, ra.xs[[2]])
plot(WSE~w, ra.xs[[3]])



