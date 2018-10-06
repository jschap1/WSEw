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
# Revised 9/27/2018
#   Consolidated/cleaned up code so there are fewer files, more abstraction
# Revised 10/4/2018
#   Fixed problems related to negative A0 values, which should only occur when
#   there is a large amount of error, so the fits slope the wrong way.
#   Removed most of the plotting code, relegated to "scrap." Will code up in make_plots later.

# ------------------------------------------------------------------------------------------------
# Set up environment

library(foreach)
library(doMC)
library(raster)
library(strucchange)
library(minpack.lm)
library(WSEw)

setwd("/Users/jschap/Desktop/Cross_Sections")
truth_dir <- "/Volumes/HD3/Cross_Sections/true_parameters"
opar <- par()

library(devtools)
library(roxygen2)
update_WSEw()

# Experiment description
nr <- 3
reach_avg <- "10km"
spacing <- 5 # m
swot_sampling <- "even"
pool <- 21
err_type <- "MC"
M <- 500 # number of replicates

expo <- seq(0.05, 0.95, length.out = 19) # exposure levels
n_exp_levels <- length(expo)
nr <- length(rWSEw)

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

set.seed(704753262) # this doesn't really set the seed properly when using foreach. need a different method.

# ------------------------------------------------------------------------------------------------
# Load data from previous runs (optional)

# Cross section and WSE-w data
load(file.path(truth_dir, "processed_xs_data_10km.rda"))
# loads {cross_sections, cross_sections_avg, xWSEw, rWSEw}

# Fitted models
lf <- readRDS(file.path(exp_dir, "lf.rds"))
sb <- readRDS(file.path(exp_dir, "sb.rds")) 

# True hydraulic parameters
# w0.ra <- readRDS(file.path(truth_dir, "w0_ra.rds"))
# WP.true.xs <- readRDS(file.path(truth_dir, "WP_true_xs.rds"))
WP.true.ra <- readRDS(file.path(truth_dir, "WP_true_ra.rds")) 

# A.true.xs <- readRDS(file.path(truth_dir, "A_true_xs.rds"))
A.true.ra <- readRDS(file.path(truth_dir, "A_true_ra.rds"))
A0.true.ra <- readRDS(file.path(truth_dir, "A0_true_ra.rds"))

s0.true.ra <- readRDS(file.path(truth_dir, "s0_true_ra.rds"))
# s0.true.xs <- readRDS(file.path(truth_dir, "s0_true_xs.rds"))
z0.true.ra <- readRDS(file.path(truth_dir, "z0_true_ra.rds"))
# z0.true.xs <- readRDS(file.path(truth_dir, "z0_true_xs.rds"))

# Predicted hydraulic parameters
load(file.path(exp_dir, "pred_lf_tmp.rda"))
load(file.path(exp_dir, "pred_sb_tmp.rda"))
load(file.path(exp_dir, "pred_nl_tmp.rda"))
load(file.path(exp_dir, "pred_nlsb_tmp.rda"))

load(file.path(exp_dir, "z0_pred.rda"))
load(file.path(exp_dir, "A_pred.rda"))
load(file.path(exp_dir, "WP_pred.rda"))
load(file.path(exp_dir, "A0_pred.rda"))
load(file.path(exp_dir, "s0_pred.rda"))

# Prediction statistics
# load(file.path(exp_dir, "bias.rda"))
# load(file.path(exp_dir, "rmse.rda"))

# ------------------------------------------------------------------------------------------------
# UMESC data processing

refWSE <- 667 # refWSE for pool 4 (ft)
refWSE <- refWSE/(39.37/12) #  convert to m

# Bathymetry
if (!file.exists("Data/p21_depth.tif"))
{
  bathy.dir <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC"
  bathy.name <- "bath_pool_4/bath_1992_p4/w001001.adf"
  bathyfile <- file.path(bathy.dir, bathy.name)
  #levels(umesc)
  umesc <- raster(bathyfile)
  depth_5 <- as_numeric_raster(umesc, att = "DEPTH_M") # native resolution depths (5 m), takes a while (like 30 minutes or something)
  
  depth_5_refined <- depth_5 # there are some very large "depth" values; this gets rid of them
  depth_5_refined[values(depth_5>100)] <- NA
  depth_5 <- depth_5_refined
  
  writeRaster(depth_5, file = "Data/p4_depth.tif", overwrite = TRUE, progress = "text")
} else 
{
  depth_5 <- raster("Data/p21_depth.tif")
  # depth_50 <- aggregate(depth_5, fact = 10) # resample depth to 50 m resolution
}

# River centerline
riv.dir <- "/Users/jschap/Desktop/Cross_Sections/Data/Centerlines"
load(file.path(riv.dir, "centerline4.rda"))
riv <- centerline_p4

# Load USGS stream gauges in UMRB
gauges <- read.table("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/Stage/gauges_UMB_QC.txt")
names(gauges) <- c("id","lat","lon","V4")
coordinates(gauges) <- ~lon + lat # make SpatialPoints object
crs(gauges) <- "+init=epsg:4326"
gauges.utm <- spTransform(gauges, crs(riv))

# Plot the study area
plot(depth_5, main = "UMRB Pool 4", xlab = "easting", ylab = "northing")
lines(riv)
points(gauges.utm, pch = 19, col = "black", cex = 0.5)

# Compute cross section data from raw bathymetry (transects_name is defunct)
cross_sections <- auto_transects(section_length = 5, depth = depth_5, refWSE = refWSE,
                                 savename = transects_name, makeplot = FALSE, riv = riv)
# (Takes about 4 hours at the highest possible resolution)

# Method 1: find corresponding flow width for WSE ranging from empty to bankfull conditions
# xWSEw <- calc_WSEw(cross_sections, interval = 0.05, dx = 1) # number of data points depends on discretization

# rWSEw <- reach_avg(xWSEw)
# rWSEw_burr <- reach_avg(xWSEw)

# save(cross_sections, xWSEw, rWSEw, file = file.path(truth_dir, "processed_xs_data_10km.rda"))

# Plot the observations
# plot(WSE~w, rWSEw[[1]], main = "WSE-w sampling, three years")
# points(WSE~w, rWSEw_burr[[1]], col="red", pch=19)
# legend("topleft", legend = c("Even sampling","Burr sampling"), fill = c("black", "red"))

# Make reach-average effective cross sections that are 10 km long each
# hist(unlist(lapply(cross_sections$x, length)))
xs.res <- resample_xs(cross_sections, n = 300)
xs.avg <- vector(length = 3, "list")
xs.avg[[1]] <- calc_mean_cross_section(xs.res[1:2000,,])
xs.avg[[2]] <- calc_mean_cross_section(xs.res[2001:4000,,])
xs.avg[[3]] <- calc_mean_cross_section(xs.res[4001:5773,,])

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

save(cross_sections, cross_sections_avg, xWSEw, rWSEw, file = file.path(truth_dir, "processed_xs_data_10km.rda"))

par(mfrow = c(1,2))
plot(b~x, xs.avg[[1]], type = "l", 
     main = "Reach average cross section 1",
     xlim = c(0,800), ylim = c(135,144))
plot(WSE~w, rWSEw[[1]], xlab = "width (m)", ylab = "height (m)", 
     ylim = c(135, 144), main = "Height-width relationship", type = "o")
plot(b~x, xs.avg[[2]], type = "l", 
     main = "Reach average cross section 2",
     xlim = c(0,650), ylim = c(135,144))
plot(WSE~w, rWSEw[[2]], xlab = "width (m)", ylab = "height (m)", 
     ylim = c(135, 144), main = "Height-width relationship", type = "o")
plot(b~x, xs.avg[[3]], type = "l", 
     main = "Reach average cross section 3",
     xlim = c(0,800), ylim = c(135,144))
plot(WSE~w, rWSEw[[3]], xlab = "width (m)", ylab = "height (m)", 
     ylim = c(135, 144), main = "Height-width relationship", type = "o")

# ------------------------------------------------------------------------------------------------
# Calculate true parameters

# z0
z0.true.xs <- unlist(lapply(cross_sections$b, min))
z0.true.ra <- unlist(lapply(cross_sections_avg$b, min))
saveRDS(z0.true.xs, file = file.path(truth_dir, "z0_true_xs.rds"))
saveRDS(z0.true.ra, file = file.path(truth_dir, "z0_true_ra.rds")) # cross section

# s0
s0.true.xs <- diff(z0.true.xs)
s0.true.ra <- diff(z0.true.ra)
par(mfrow=c(2,1))
plot(s0.true.xs, type = "l")
plot(s0.true.ra, type = "l")
abline(0,0)
saveRDS(s0.true.xs, file = file.path(truth_dir, "s0_true_xs.rds"))
saveRDS(s0.true.ra, file = file.path(truth_dir, "s0_true_ra.rds"))

# A, WP
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
saveRDS(A.true.ra, file = file.path(truth_dir, "A_true_ra.rds"))
saveRDS(WP.true.ra, file = file.path(truth_dir, "WP_true_ra.rds"))

# A0
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

# A0 (reach average)
A0.true.ra <- array(dim = c(nr, n_exp_levels))
# w1.ra <- array(dim = c(nr, n_exp_levels)) # lowest observed width value
# h1.ra <- array(dim = c(nr, n_exp_levels)) # lowest observed height value
for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    WSEw_obs <- observe(rWSEw[[r]], sd_wse = 0, sd_w = 0, exposure = expo[k])
    p <- max(which(rWSEw[[r]]$WSE<min(WSEw_obs$WSE))) # index up to which is not observed
    # w0.ra[r,k] <- min(WSEw_obs$w)
    # h1.ra[r,k] <- min(WSEw_obs$WSE)
    A0.true.ra[r,k] <- calc_A_from_WSEw(rWSEw[[r]][1:(p+1),])
  }
  print(r)
}
# saveRDS(w0.ra, file = file.path(exp_dir, "w0_ra.rds"))
# saveRDS(h1.ra, file = file.path(exp_dir, "h1_ra.rds"))
saveRDS(A0.true.ra, file = file.path(truth_dir, "A0_true_ra.rds"))

# ------------------------------------------------------------------------------------------------
# Fit WSE-w models

# Model fitting at different exposure levels, known measurements
# Should use errors in variables method, but this can be refined later.
# Will do this with uncertain measurements later, preferably without MC simulation

# Run the parallel computations (see computation_time_calculator.xls)

registerDoMC(cores = 3)
# ncores <- detectCores()
# registerDoMC(cores = ncores - 1)

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

# ------------------------------------------------------------------------------------------------
# Predict hydraulic parameters

# Do computations in parallel
# ncores <- detectCores()
# registerDoMC(cores = ncores - 1)
registerDoMC(cores = 3)

begin.time <- Sys.time()
pred_lf <- foreach(r = 1:nr) %dopar% {pred_linear_par(r)}
print(Sys.time() - begin.time)
save(pred_lf, file = file.path(exp_dir, "pred_lf_tmp.rda"))

begin.time <- Sys.time()
pred_sb <- foreach(r = 1:nr) %dopar% {pred_sb_par(r)}
print(Sys.time() - begin.time)
save(pred_sb, file = file.path(exp_dir, "pred_sb_tmp.rda"))

# for speed, it is best to pre-load w1, h1
# This is done below
h1 <- array(dim = c(nr, n_exp_levels, M))
w1 <- array(dim = c(nr, n_exp_levels, M))
for (r in 1:nr)
{
  obsname <- paste0("obs/WSEw_obs_r_", r, ".rds") # load observations for this reach
  WSEw_obs <- readRDS(file.path(exp_dir, obsname))
  for (m in 1:M)
  {
    for (k in 1:n_exp_levels)
    {
      h1[r,k,m] <- WSEw_obs[[k]][[m]]$WSE[1]
      w1[r,k,m] <- WSEw_obs[[k]][[m]]$w[1]
    }
  }
}
saveRDS(h1, file.path(exp_dir, "h1_obs.rds"))
saveRDS(w1, file.path(exp_dir, "w1_obs.rds"))

begin.time <- Sys.time()
pred_nl <- foreach(r = 1:nr) %dopar% {pred_nl_par(r, rWSEw, w1 = w1, h1 = h1)}
print(Sys.time() - begin.time)
save(pred_nl, file = file.path(exp_dir, "pred_nl_tmp.rda"))

begin.time <- Sys.time()
pred_nlsb <- foreach(r = 1:nr) %dopar% {pred_nlsb_par(r, rWSEw, w1 = w1, h1 = h1)}
print(Sys.time() - begin.time)
save(pred_nlsb, file = file.path(exp_dir, "pred_nlsb_tmp.rda"))

# Sometimes there are negative A0 values.
# This occurs when the z0 prediction is higher than the observed h1 and can be attributed to measurement error.
# This is physically inconsistent, so neglect these predictions and analyze the remaining predicted values.

pred_lf <- clean_by_A0(pred_lf, h1) # remove A0 and z0 predictions where A0<0
pred_sb <- clean_by_A0(pred_sb, h1)
pred_nl <- clean_by_A0(pred_nl, h1)
pred_nlsb <- clean_by_A0(pred_nlsb, h1)

# Run predictions again for any reaches where errors occurred

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

# !!!
# names(z0.l) # coerce to data frames and save as sample data for the R package.

# Compute slope via finite difference
s0.l <- apply(z0.l, c(2,3), diff)
s0.sb <- apply(z0.sb, c(2,3), diff)
s0.nl <- apply(z0.nl, c(2,3), diff)
s0.nlsb <- apply(z0.nlsb, c(2,3), diff)
save(s0.l, s0.sb, s0.nl, s0.nlsb, file = file.path(exp_dir, "s0_pred.rda"))

# ----------------------------------------------------------------------------------------------------
# Calculate prediction error

z0.l.error <- calc_prediction_error(z0.l, z0.true.ra)
z0.sb.error <- calc_prediction_error(z0.sb, z0.true.ra)
z0.nl.error <- calc_prediction_error(z0.nl, z0.true.ra)
z0.nlsb.error <- calc_prediction_error(z0.nlsb, z0.true.ra)

A0.l.error <- calc_prediction_error(A0.l, A0.true.ra, TRUE)
A0.sb.error <- calc_prediction_error(A0.sb, A0.true.ra, TRUE)
A0.nl.error <- calc_prediction_error(A0.nl, A0.true.ra, TRUE)
A0.nlsb.error <- calc_prediction_error(A0.nlsb, A0.true.ra, TRUE)
