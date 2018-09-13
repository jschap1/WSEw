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

# ------------------------------------------------------------------------------------------------

# Set up environment
library(foreach)
library(doMC)
library(raster)
library(strucchange)
library(segmented)
library(minpack.lm)
library(actuar) # needed if using Burr distribution in calc_WSEw3
library(WSEw)

setwd("/Users/jschap/Desktop/Cross_Sections")
registerDoMC(cores = 7)

# Experiment description
n_exp_levels <- 19
nr <- 50
reach_avg <- TRUE
spacing <- 5 # m
swot_sampling <- "even"
pool <- 21
err_type <- "mc"
M <- 100 # number of replicates

# Make a directory to store results
exp_desc <- paste0("pool_", pool, "_ra_",reach_avg,"_nr_",nr,"_expo_",n_exp_levels,"_spacing_",spacing,"_sampling_",swot_sampling, "_", err_type, "_replicates_", M)
fits_dir <- "/Users/jschap/Desktop/Cross_Sections/Outputs/" # directory for modeling outputs
exp_dir <- file.path(fits_dir, exp_desc) # directory for this experiment's outputs

if (!dir.exists(exp_dir))
{
  dir.create(exp_dir)
  dir.create(file.path(exp_dir, "lf"))
  dir.create(file.path(exp_dir, "sb"))
  dir.create(file.path(exp_dir, "sbm"))
  dir.create(file.path(exp_dir, "nl"))
  dir.create(file.path(exp_dir, "nlsb"))
}

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
sbm <- readRDS(file.path(exp_dir, "sbm.rds"))

nl1 <- readRDS(file.path(exp_dir, "nl_1_to_3500.rds")) # ~ 2 min, 8.8 Gb
nl2 <- readRDS(file.path(exp_dir, "nl3501to3774.rds")) # 16 seconds, 0.6 Gb
nl <- c(nl1, nl2)
rm(nl1, nl2)
gc()

begin.time <- Sys.time() # 6.5 min, ~30 Gb
nlsb1 <- readRDS(file.path(exp_dir, "nlsb1to500.rds"))
nlsb2 <- readRDS(file.path(exp_dir, "nlsb501to1000.rds"))
nlsb3 <- readRDS(file.path(exp_dir, "nlsb1001to1500.rds"))
nlsb4 <- readRDS(file.path(exp_dir, "nlsb1501to2000.rds"))
nlsb5 <- readRDS(file.path(exp_dir, "nlsb2001to2500.rds"))
nlsb6 <- readRDS(file.path(exp_dir, "nlsb2501to3000.rds"))
nlsb7 <- readRDS(file.path(exp_dir, "nlsb3001to3500.rds"))
nlsb8 <- readRDS(file.path(exp_dir, "nlsb3501to3774.rds"))
duration <- Sys.time() - begin.time
print(duration)

nlsb <- c(nlsb1, nlsb2, nlsb3, nlsb4, nlsb5, nlsb6, nlsb7, nlsb8) # this is 19.1 Gb
rm(nlsb1, nlsb2, nlsb3, nlsb4, nlsb5, nlsb6, nlsb7, nlsb8)
gc()

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
WSEw_obs <- readRDS(file.path(exp_dir, "WSEw_obs.rds"))

# ------------------------------------------------------------------------------------------------
# Load predicted hydraulic parameters

load(file.path(exp_dir, "z0_pred.rda"))
load(file.path(exp_dir, "A_pred.rda"))
load(file.path(exp_dir, "WP_pred.rda"))
load(file.path(exp_dir, "A0_pred.rda"))

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

# Method 3: use a probability distribution based on stage data to simulate SWOT observations
xWSEw1 <- calc_WSEw3(cross_sections, dist = "burr", n.obs = floor(1*365/10)) # one year of obs
xWSEw3 <- calc_WSEw3(cross_sections, dist = "burr", n.obs = floor(1*365/10)) # 3 yrs of obs

rWSEw <- reach_avg(xWSEw1, l = 10000, res = 5)

save(cross_sections, xWSEw, rWSEw, file = file.path(saveloc, "/Data/Processed_Data/processed_xs_data.rda"))

# Plot the observations
plot(WSE~w, rWSEw[[3000]], main = "WSE-w sampling, one year")
points(WSE~w, xWSEw1[[1000]], col="red", pch=19)
legend("topleft", legend = c("Even sampling","Burr sampling"), fill = c("black", "red"))

# ------------------------------------------------------------------------------------------------
# Main experiments - model fitting at different exposure levels, known measurements
# Should use errors in variables method, but this can be refined later.
# Will do this with uncertain measurements later, preferably without MC simulation

expo <- seq(0.05, 0.95, length.out = n_exp_levels) # exposure levels
n_exp_levels <- length(expo)
nr <- length(rWSEw)

# Initialize WSEw_obs
WSEw_obs <- vector(length = nr, "list")
for (r in 1:nr)
{
  WSEw_obs[[r]] <- vector(length = n_exp_levels, "list")
}
for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    WSEw_obs[[r]][[k]] <- vector(length = M, "list")
  }
}

# ---------------------------------------------------------------------------------------------------------------

# Make observations
set.seed(704753262)
begin.time <- Sys.time() # takes 41 s for nr = 50, M = 100
for (r in 1:nr) # loop over reaches
{
  for (k in 1:n_exp_levels) # loop over exposure levels
  {
    for (m in 1:M)
    {
      WSEw_obs[[r]][[k]][[m]] <- observe(WSEw = rWSEw[[r]], exposure = expo[k])
    }
    print(k)
  }
  ## save data for cross section (to avoid bulky files)
  #fname <- paste0("WSEw_obs_r_", r, ".rds")
  #saveRDS(WSEw_obs, file = file.path(exp_dir, fname))
  if (r%%5 == 0)
  {
    current.time <- Sys.time()
    te <- 
      print(paste("Processed", r, "of", nr, "cross sections")) # display progress
    print(current.time - begin.time)
  }
}
saveRDS(WSEw_obs, file = file.path(exp_dir, "WSEw_obs.rds"))

source("Codes/modelfit_par.R") # load the fitting functions

# Run the parallel computations (see computation_time_calculator.xls)
lfval <- foreach(r = 1:nr, .combine = c) %dopar% {fit_linear_par(r)} # returns a dummy value; the point is to save lf.rds
sbval <- foreach(r = 1:nr, .combine = c) %dopar% {fit_sb_par(r)}
sbmval <- foreach(r = 1:nr, .combine = c) %dopar% {fit_sbm_par(r)}
nlval <- foreach(r = 1:nr, .combine = c) %dopar% {fit_nl_par(r)}
nlsbval <- foreach(r = 1:nr, .combine = c) %dopar% {fit_nlsb_par(r)}

# begin.time <- Sys.time()
# nlsbval <- foreach(r = 1:nr, .combine = c) %dopar% {fit_nlsb_par(r)}
# print(Sys.time() - begin.time)

# ------------------------------------------------------------------------------------------------
# Save fitted models

# This is done as the models are fitted, above.
# Still need to revise code below for Monte Carlo implementation.

# ------------------------------------------------------------------------------------------------------------
# Compute true z0

z0.true.xs <- unlist(lapply(cross_sections$b, min))
z0.true.ra <- ra(z0.true.xs, n = 2000)
saveRDS(z0.true.xs, file = file.path(exp_dir, "z0_true_xs.rds")) # cross section
saveRDS(z0.true.ra, file = file.path(exp_dir, "z0_true_ra.rds")) # reach-averaged

# ------------------------------------------------------------------------------------------------------------
# Compute true s0

s0.true.xs <- diff(z0.true.xs)
s0.true.ra <- diff(z0.true.ra)
par(mfrow=c(2,1))
plot(s0.true.xs, type = "l")
plot(s0.true.ra, type = "l")
abline(0,0)
saveRDS(s0.true.xs, file = file.path(exp_dir, "s0_true_xs.rds")) 
saveRDS(s0.true.ra, file = file.path(exp_dir, "s0_true_ra.rds"))

# ------------------------------------------------------------------------------------------------------------
# Compute true A, WP

# Calculate true hydraulic parameters for cross sections
n.xs <- length(xWSEw)
A.true.xs <- vector(length = n.xs)
WP.true.xs <- vector(length = n.xs)
for (k in 1:n.xs)
{
  x <- cross_sections$x[[k]]
  b <- cross_sections$b[[k]]
  WSE <- max(b) # bankfull
  A.true.xs[k] <- calc_A(x, b, WSE)
  WP.true.xs[k] <- calc_WP(x, b, WSE)
}
plot(A.true.xs, type="l")
plot(WP.true.xs, type="l")
saveRDS(A.true.xs, file = file.path(exp_dir, "A_true_xs.rds")) 
saveRDS(WP.true.xs, file = file.path(exp_dir, "WP_true_xs.rds")) 

# Calculate reach-average hydraulic parameters
n <- 2000 # number of segments in a reach = 10 km/5 m = 2000
A.true.ra <- ra(A.true.xs, n)
WP.true.ra <- ra(WP.true.xs, n)
plot(A.true.ra, type="l")
plot(WP.true.ra, type="l")
saveRDS(A.true.ra, file = file.path(exp_dir, "A_true_ra.rds")) 
saveRDS(WP.true.ra, file = file.path(exp_dir, "WP_true_ra.rds")) 

# ------------------------------------------------------------------------------------------------------------
# Compute true A0

# Calculate true A0
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

# Initialize
z0.l <- array(dim = c(nr, n_exp_levels, M))
z0.sb <- array(dim = c(nr, n_exp_levels, M))
z0.sbm <- array(dim = c(nr, n_exp_levels, M))
z0.nl <- array(dim = c(nr, n_exp_levels, M))
z0.nlsb <- array(dim = c(nr, n_exp_levels, M))

A.l <- array(dim = c(nr, n_exp_levels, M))
A.sb <- array(dim = c(nr, n_exp_levels, M))
A.sbm <- array(dim = c(nr, n_exp_levels, M))
A.nl <- array(dim = c(nr, n_exp_levels, M))
A.nlsb <- array(dim = c(nr, n_exp_levels, M))

WP.l <- array(dim = c(nr, n_exp_levels, M))
WP.sb <- array(dim = c(nr, n_exp_levels, M))
WP.sbm <- array(dim = c(nr, n_exp_levels, M))
WP.nl <- array(dim = c(nr, n_exp_levels, M))
WP.nlsb <- array(dim = c(nr, n_exp_levels, M))

A0.l <- array(dim = c(nr, n_exp_levels, M))
A0.sb <- array(dim = c(nr, n_exp_levels, M))
A0.sbm <- array(dim = c(nr, n_exp_levels, M))
A0.nl <- array(dim = c(nr, n_exp_levels, M))
A0.nlsb <- array(dim = c(nr, n_exp_levels, M))

# Linear
for (r in 1:nr) # takes 5 minutes for 50 cross sections at 19 exposure levels
{
  lf_name <- paste0("lf/lf_", "r_", r, ".rds")
  lf <- readRDS(file.path(exp_dir, lf_name))
  for (k in 1:n_exp_levels)
  { 
    for (m in 1:M)
    {
      # if statements handle cases where the model is NULL/no model was fit
      if (class(lf[[k]][[m]]) == "lm")
      {
        z0.l[r,k,m] <- predict(lf[[k]][[m]], newdata = data.frame(w = 0))
        A.l[r,k,m] <- calc_model_A(lf[[k]][[m]], type = "linear")
        WP.l[r,k,m] <- calc_model_WP(lf[[k]][[m]], type = "linear")
        A0.l[r,k,m] <- calc_model_A0(lf[[k]][[m]], type = "linear")
      }
    }
  }
  if (r %% 1 == 0)
  {
    print(paste("progress:", r, "of", nr))
  }
}

# Slope break
for (r in 1:nr)
{
  sb_name <- paste0("sb/sb_", "r_", r, ".rds")
  sb <- readRDS(file.path(exp_dir, sb_name))
  for (k in 1:n_exp_levels)
  { 
    for (m in 1:M)
    {
      if (class(sb[[k]][[m]][[1]]) == "lm")
      {
        z0.sb[r,k,m] <- predict(sb[[k]][[m]][[1]], newdata = data.frame(w = 0))
        A.sb[r,k,m] <- calc_model_A(sb[[k]][[m]], type = "sb")
        WP.sb[r,k,m] <- calc_model_WP(sb[[k]][[m]], type = "sb")
        A0.sb[r,k,m] <- calc_model_A0(sb[[k]][[m]], type = "sb")
      }
    }
  }
  if (r %% 1 == 0)
  {
    print(paste("progress:", r, "of", nr))
  }
}

# SBM
for (r in 1:nr)
{
  sbm_name <- paste0("sbm/sbm_", "r_", r, ".rds")
  sbm <- readRDS(file.path(exp_dir, sbm_name))
  for (k in 1:n_exp_levels)
  { 
    for (m in 1:M)
    {
      if (any(class(sbm[[k]][[m]][[1]])=="lm"))
      {
        z0.sbm[r,k,m] <- predict(sbm[[k]][[m]][[1]], newdata = data.frame(w = 0))
        A.sbm[r,k,m] <- calc_model_A(sbm[[k]][[m]], type = "sbm")
        WP.sbm[r,k,m] <- calc_model_WP(sbm[[k]][[m]], type = "sbm")
        A0.sbm[r,k,m] <- calc_model_A0(sbm[[k]][[m]], type = "sbm")
      }
    }
  }
  if (r %% 1 == 0)
  {
    print(paste("progress:", r, "of", nr))
  }
}

# Nonlinear
for (r in 1:nr)
{
  nl_name <- paste0("nl/nl_", "r_", r, ".rds")
  nl <- readRDS(file.path(exp_dir, nl_name))
  for (k in 1:n_exp_levels)
  { 
    for (m in 1:M)
    {
      if (class(nl[[k]][[m]]) == "nls")
      {
        z0.nl[r,k,m] <- predict(nl[[k]][[m]], newdata = data.frame(w = 0))
        A.nl[r,k,m] <- calc_model_A(nl[[k]][[m]], type = "nl", WSEw = rWSEw[[r]])
        WP.nl[r,k,m] <- calc_model_WP(nl[[k]][[m]], type = "nl", w = rWSEw[[r]]$w)
        A0.nl[r,k,m] <- calc_model_A0(nl[[k]][[m]], type = "nl", w0 = w0.ra[r,k])
      }
    }
  }
  if (r %% 1 == 0)
  {
    print(paste("progress:", r, "of", nr))
  }
}

# NLSB
for (r in 1:nr)
{
  nlsb_name <- paste0("nlsb/nlsb_", "r_", r, ".rds")
  nlsb <- readRDS(file.path(exp_dir, nlsb_name))
  for (k in 1:n_exp_levels)
  { 
    for (m in 1:M)
    {
      if (class(nlsb[[k]][[m]][[1]]) == "nls")
      {
        z0.nlsb[r,k,m] <- predict(nlsb[[k]][[m]][[1]], newdata = data.frame(w = 0))
        A.nlsb[r,k,m] <- calc_model_A(nlsb[[k]][[m]], type = "nlsb", WSEw = rWSEw[[r]]) # there may be a bug in the type = nlsb code here
        WP.nlsb[r,k,m] <- calc_model_WP(nlsb[[k]][[m]], type = "nlsb", w = rWSEw[[r]]$w)
        A0.nlsb[r,k,m] <- calc_model_A0(nlsb[[k]][[m]], type = "nlsb", w0 = w0.ra[r,k])
      }
    }
  }
  if (r %% 1 == 0)
  {
    print(paste("progress:", r, "of", nr))
  }
}

save(z0.l, z0.sb, z0.sbm, z0.nl, z0.nlsb, file = file.path(exp_dir, "z0_pred.rda"))
save(A.l, A.sb, A.sbm, A.nl, A.nlsb, file = file.path(exp_dir, "A_pred.rda"))
save(WP.l, WP.sb, WP.sbm, WP.nl, WP.nlsb, file = file.path(exp_dir, "WP_pred.rda"))
save(A0.l, A0.sb, A0.sbm, A0.nl, A0.nlsb, file = file.path(exp_dir, "A0_pred.rda"))

# Compute slope via finite difference
s0.l <- apply(z0.l, c(2,3), diff)
s0.sb <- apply(z0.sb, c(2,3), diff)
s0.sbm <- apply(z0.sbm, c(2,3), diff)
s0.nl <- apply(z0.nl, c(2,3), diff)
s0.nlsb <- apply(z0.nlsb, c(2,3), diff)
save(s0.l, s0.sb, s0.sbm, s0.nl, s0.nlsb, file = file.path(exp_dir, "s0_pred.rda"))

# ------------------------------------------------------------------------------------------------
# Calculate A0 prediction error for linear model with no measurement error

# A0.var <- array(dim = c(nr, n_exp_levels))
# for (r in 1:nr)
# {
#   for (k in 1:n_exp_levels)
#   {
#     if (class(lf[[r]][[k]]) == "lm")
#     {
#       A0.var[r,k] <- calc_A0_variance(lf[[r]][[k]])
#     }
#   }
# }
# A0.sd <- sqrt(A0.var)
# summary(A0.sd)
# saveRDS(A0.sd, file.path(exp_dir, "A0_sd.rds"))

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
# This is as far as I've gotten.
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

z0.sbm.lower <- array(dim = c(nr, n_exp_levels))
z0.sbm.med <- array(dim = c(nr, n_exp_levels))
z0.sbm.upper <- array(dim = c(nr, n_exp_levels))

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
    
    q1 <- quantile(z0.sbm[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    z0.sbm.lower[r,k] <- as.numeric(q1[1])
    z0.sbm.med[r,k] <- as.numeric(q1[2])
    z0.sbm.upper[r,k] <- as.numeric(q1[3])
    
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

A0.sbm.lower <- array(dim = c(nr, n_exp_levels))
A0.sbm.med <- array(dim = c(nr, n_exp_levels))
A0.sbm.upper <- array(dim = c(nr, n_exp_levels))

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
    
    q1 <- quantile(A0.sbm[r,k,], c(0.05, 0.5, 0.95), na.rm = TRUE)
    A0.sbm.lower[r,k] <- as.numeric(q1[1])
    A0.sbm.med[r,k] <- as.numeric(q1[2])
    A0.sbm.upper[r,k] <- as.numeric(q1[3])
    
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

# ------------------------------------------------------------------------------------------------
# Make a filled plot showing percentiles of the A0 predictions for each estimate

# base r example:
age <- 1:10
y.low <- rnorm(length(age), 150, 25) + 10*age
y.high <- rnorm(length(age), 250, 25) + 10*age

plot(age,y.high,type = 'n', ylim = c(100, 400),
     ylab = 'Y Range', xlab = 'Age (years)')
lines(age, y.low, col = 'grey')
lines(age, y.high, col = 'grey')

polygon(c(age, rev(age)), c(y.high, rev(y.low)),
        col = "grey30", border = NA)

# going to have better luck with ggplot2. Learn it.


# ------------------------------------------------------------------------------------------------
# Make plots of z0, A, WP, A0, s0 along the river

k <- 8 # exposure level
plot(z0.true.ra, main = paste("z0 (m) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(120,138), xlab = "reach-average cross section", ylab = "min. bed elevation")
lines(z0.l[,k], col = "red")
lines(z0.sb[,k], col = "orange")
lines(z0.sbm[,k], col = "purple")
lines(z0.nl[,k], col = "green")
lines(z0.nlsb[,k], col = "blue")
legend("bottomleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

plot(A.true.ra, main = paste("Flow area (m^2) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(0,6000))
lines(A.l[,k], col = "red")
lines(A.sb[,k], col = "orange")
lines(A.sbm[,k], col = "purple")
lines(A.nl[,k], col = "green")
lines(A.nlsb[,k], col = "blue")
legend("bottomleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

plot(WP.true.ra, main = paste("Wetted perimeter (m) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(450,850))
lines(WP.l[,k], col = "red")
lines(WP.sb[,k], col = "orange")
lines(WP.sbm[,k], col = "purple")
lines(WP.nl[,k], col = "green")
lines(WP.nlsb[,k], col = "blue")
legend("topleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# A0
plot(A0.true.ra[,k], main = paste("A0 (m^2) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(0,4500))
lines(A0.l[,k], col = "red")
lines(A0.sb[,k], col = "orange")
lines(A0.sbm[,k], col = "purple")
lines(A0.nl[,k], col = "green")
lines(A0.nlsb[,k], col = "blue")
legend("topleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# s0
plot(s0.true.ra*1e5, main = paste("s0 (cm/km) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(-2000,2000), 
     xlab = "reach-average cross section", ylab = "bed slope")
lines(1e5*s0.l[,k], col = "red")
lines(1e5*s0.sb[,k], col = "orange")
lines(1e5*s0.sbm[,k], col = "purple")
lines(1e5*s0.nl[,k], col = "green")
lines(1e5*s0.nlsb[,k], col = "blue")
legend("topleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

K <- 1000 # Smooth slope using a k-point moving average
s0.true.ra.smooth <- filter(s0.true.ra, sides = 2, filter = rep(1/K,K))
s0.l.smooth <- filter(s0.l, sides = 2, filter = rep(1/K,K))
s0.sb.smooth <- filter(s0.sb, sides = 2, filter = rep(1/K,K))
s0.sbm.smooth <- filter(s0.sbm, sides = 2, filter = rep(1/K,K))
s0.nl.smooth <- filter(s0.nl, sides = 2, filter = rep(1/K,K))
s0.nlsb.smooth <- filter(s0.nlsb, sides = 2, filter = rep(1/K,K))

# s0 smoothed
plot(s0.true.ra.smooth*1e5, main = paste("s0 (cm/km) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(-500,500), 
     xlab = "reach-average cross section", ylab = "smoothed bed slope")
lines(1e5*s0.l.smooth[,k], col = "red")
lines(1e5*s0.sb.smooth[,k], col = "orange")
lines(1e5*s0.sbm.smooth[,k], col = "purple")
lines(1e5*s0.nl.smooth[,k], col = "green")
lines(1e5*s0.nlsb.smooth[,k], col = "blue")
legend("topleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# ------------------------------------------------------------------------------------------------
# Make plots of average z0, A, WP, s0, A0 error at each exposure level

# z0
pred_z0 <- list(z0.l, z0.sb, z0.sbm, z0.nl, z0.nlsb)
z0.bias <- plot_bias(expo, pred_z0, z0.true, na.rm = TRUE,
                     main = "z0 bias vs. exposure level, no meas. error", ylab = "Bias (m)", ylim = c(-25,3))
legend("bottomright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# s0
pred.s0.cmkm <- list(1e5*s0.l, 1e5*s0.sb, 1e5*s0.sbm, 1e5*s0.nl, 1e5*s0.nlsb)
s0.true.cmkm <- s0.true.ra*1e5 # convert units to cm/km
s0.bias <- plot_bias(expo, pred.s0.cmkm, s0.true.cmkm, na.rm = TRUE,
                     main = "s0 bias vs. exposure level, no meas. error", ylab = "Bias (cm/km)", ylim = c(-50,400))
legend("topright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# A
pred_A <- list(A.l, A.sb, A.sbm, A.nl, A.nlsb)
A.bias <- plot_bias(expo, pred_A, A.true.ra, na.rm = TRUE,
                    main = "A bias vs. exposure level, no meas. error", ylab = "Bias (m^2)")
legend("topright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# WP
pred_WP <- list(WP.l, WP.sb, WP.sbm, WP.nl, WP.nlsb)
WP.bias <- plot_bias(expo, pred_WP, WP.true.ra, na.rm = TRUE,
                     main = "WP bias vs. exposure level, no meas. error", ylab = "Bias (m)")
legend("topright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# A0
A0.bias <- array(dim = c(n_exp_levels, 5)) # cannot use plot_bias because A0 changes with exposure level
na.rm = TRUE
for (k in 1:n_exp_levels)
{
  A0.bias[k,1] <- mean(A0.l[,k] - A0.true[,k], na.rm = na.rm)
  A0.bias[k,2] <- mean(A0.sb[,k] - A0.true[,k], na.rm = na.rm)
  A0.bias[k,3] <- mean(A0.sbm[,k] - A0.true[,k], na.rm = na.rm)
  A0.bias[k,4] <- mean(A0.nl[,k] - A0.true[,k], na.rm = na.rm)
  A0.bias[k,5] <- mean(A0.nlsb[,k] - A0.true[,k], na.rm = na.rm)
}
A0.bias <- as.data.frame(A0.bias)
names(A0.bias) <- c("l","sb","sbm","nl","nlsb")

# Plot A0.bias vs. exposure level
plot(100*expo, A0.bias$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "A0.bias (m2), no meas. error", ylab = "A0.bias (m)", xlim = c(40,100), ylim = c(0,1000))
lines(100*expo, A0.bias$sb, col = "orange")
lines(100*expo, A0.bias$sbm, col = "purple")
lines(100*expo, A0.bias$nl, col = "green")
lines(100*expo, A0.bias$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

save(z0.bias, s0.bias, A.bias, WP.bias, A0.bias, file = file.path(exp_dir, "bias.rda"))

# ------------------------------------------------------------------------------------------------
# Compute RMSE

z0.rmse <- compute_rmse(pred_z0, z0.true.ra)
A.rmse <- compute_rmse(pred_A, A.true.ra)
WP.rmse <- compute_rmse(pred_WP, WP.true.ra)
s0.rmse <- compute_rmse(pred.s0.cmkm, s0.true.cmkm)

A0.rmse <- array(dim = c(n_exp_levels, 5)) # cannot use compute_rmse because A0 changes with exposure level
na.rm = TRUE
for (k in 1:n_exp_levels)
{
  A0.rmse[k,1] <- ((1/length(A0.true[,k]))*t(A0.l[,k] - A0.true[,k])%*%(A0.l[,k] - A0.true[,k]))^0.5
  A0.rmse[k,2] <- ((1/length(A0.true[,k]))*t(A0.sb[,k] - A0.true[,k])%*%(A0.sb[,k] - A0.true[,k]))^0.5
  A0.rmse[k,3] <- ((1/length(A0.true[,k]))*t(A0.sbm[,k] - A0.true[,k])%*%(A0.sbm[,k] - A0.true[,k]))^0.5
  A0.rmse[k,4] <- ((1/length(A0.true[,k]))*t(A0.nl[,k] - A0.true[,k])%*%(A0.nl[,k] - A0.true[,k]))^0.5
  A0.rmse[k,5] <- ((1/length(A0.true[,k]))*t(A0.nlsb[,k] - A0.true[,k])%*%(A0.nlsb[,k] - A0.true[,k]))^0.5
}
A0.rmse <- as.data.frame(A0.rmse)
names(A0.rmse) <- c("l","sb","sbm","nl","nlsb")

save(z0.rmse, s0.rmse, A.rmse, WP.rmse, A0.rmse, file = file.path(exp_dir, "rmse.rda"))

# ------------------------------------------------------------------------------------------------
# Make plots of z0, A, WP, s0, A0 RMSE at each exposure level

# z0
plot(100*expo, z0.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "z0 RMSE (cm/km)", ylab = "z0.rmse (m)",
     ylim = c(0,10))
lines(100*expo, z0.rmse$sb, col = "orange")
lines(100*expo, z0.rmse$sbm, col = "purple")
lines(100*expo, z0.rmse$nl, col = "green")
lines(100*expo, z0.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# s0
plot(100*expo, s0.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "s0 RMSE (cm/km)", ylab = "s0.rmse (m)",
     ylim = c(0,20000))
lines(100*expo, s0.rmse$sb, col = "orange")
lines(100*expo, s0.rmse$sbm, col = "purple")
lines(100*expo, s0.rmse$nl, col = "green")
lines(100*expo, s0.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# A
plot(100*expo, A.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "A RMSE (m2)", ylab = "A.rmse (m)",
     ylim = c(0,7000))
lines(100*expo, A.rmse$sb, col = "orange")
lines(100*expo, A.rmse$sbm, col = "purple")
lines(100*expo, A.rmse$nl, col = "green")
lines(100*expo, A.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# WP
plot(100*expo, WP.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "WP RMSE (m)", ylab = "WP.rmse (m)",
     ylim = c(0,4))
lines(100*expo, WP.rmse$sb, col = "orange")
lines(100*expo, WP.rmse$sbm, col = "purple")
lines(100*expo, WP.rmse$nl, col = "green")
lines(100*expo, WP.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# A0
plot(100*expo, A0.rmse$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "A0.rmse (m2)", ylab = "A0.rmse (m)", xlim = c(40,100), ylim = c(0,1000))
lines(100*expo, A0.rmse$sb, col = "orange")
lines(100*expo, A0.rmse$sbm, col = "purple")
lines(100*expo, A0.rmse$nl, col = "green")
lines(100*expo, A0.rmse$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)





