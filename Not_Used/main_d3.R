# Executive file for WSEw project
# 
# Created 8/2/2018 JRS
# Revised 8/10/2018 JRS
# Revised 9/5/2018 JRS 
#   Reorganized for ease of use

# ------------------------------------------------------------------------------------------------

# Set up environment
library(raster)
library(strucchange)
library(segmented)
library(minpack.lm)
library(WSEw)

setwd("Desktop/Cross_Sections")

# Experiment description
n_exp_levels <- 19
nr <- 3774
reach_avg <- TRUE
spacing <- 5 # m
swot_sampling <- "even"
pool <- 21
err_type <- "no_err"

# Make a directory to store results
exp_desc <- paste0("pool_", pool, "_ra_",reach_avg,"_nr_",nr,"_expo_",n_exp_levels,"_spacing_",spacing,"_sampling_",swot_sampling, "_", err_type)
fits_dir <- "/Users/jschap/Desktop/Cross_Sections/Outputs/" # directory for modeling outputs
exp_dir <- file.path(fits_dir, exp_desc) # directory for this experiment's outputs

if (!dir.exists(exp_dir))
{
  dir.create(exp_dir)
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
# Load fitted models

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

A0_true_ra <- readRDS(file.path(exp_dir, "lf.rds"))
s0_true_ra <- readRDS(file.path(exp_dir, "lf.rds"))

# ------------------------------------------------------------------------------------------------
# Load predicted hydraulic parameters

load("Outputs/nr3774_fitted_models_no_err/z0.rda")
load("Outputs/nr3774_fitted_models_no_err/A.rda")
load("WP.rda")

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

# Initialize outputs
lf <- vector(length = nr, "list")
sb <- vector(length = nr, "list")
sbm <- vector(length = nr, "list")
nl <- vector(length = nr, "list")
nlsb <- vector(length = nr, "list")
for (r in 1:nr)
{
  lf[[r]] <- vector(length = n_exp_levels, "list")
  sb[[r]] <- vector(length = n_exp_levels, "list")
  sbm[[r]] <- vector(length = n_exp_levels, "list")
  nl[[r]] <- vector(length = n_exp_levels, "list")
  nlsb[[r]] <- vector(length = n_exp_levels, "list")
}

# To do: use the advice here: https://stackoverflow.com/questions/12135400/errors-in-segmented-package-breakpoints-confusion
# This will likely allow sbm fits to work more often, by restarting multiple times
begin.time <- Sys.time() # It does about 17 cross sections per minute.
for (r in 1:nr) # loop over reaches
{
  for (k in 1:n_exp_levels) # loop over exposure levels
  {
    WSEw_obs <- observe(WSEw = rWSEw[[r]], exposure = expo[k], sd_wse = 0, sd_w = 0)
    lf[[r]][[k]] <- fit_linear(WSEw_obs)
    sb[[r]][[k]] <- fit_slopebreak(WSEw_obs, multiple_breaks = FALSE, continuity = TRUE)
    try(sbm[[r]][[k]] <- fit_slopebreak(WSEw_obs, multiple_breaks = TRUE, continuity = TRUE)) # sometimes this throws errors
    nl[[r]][[k]] <- fit_nonlinear(WSEw_obs)
    nlsb[[r]][[k]] <- fit_nlsb(WSEw_obs)
    if (r%%5 == 0)
    {
      current.time <- Sys.time()
      te <- current.time - begin.time
      print(paste("Processed", r, "of", nr, "cross sections")) # display progress
      print(paste("Time elapsed:", te, "units"))
    }
  }
}

# ------------------------------------------------------------------------------------------------
# Save fitted models

# It is a large amount of data. Saving may crash R.
saveRDS(lf, "Outputs/lf.rds")
saveRDS(sb, "Outputs/sb.rds") # takes about 3 minutes to save 2.3 GB
saveRDS(sbm, "Outputs/sbm.rds")

# Breaking down this save because it takes a long time
ind1 <- 0
for (j in seq(500, 3500, by = 500))
{
  ind1 <- j-500+1
  #print(ind1)
  nl1 <- nl[ind1:j]
  saveRDS(nl1, paste0("nl", ind1, "to", j, ".rds"))
  print(j)
}

# Breaking down this save because it takes a long time
ind1 <- 0
for (j in seq(500, 3500, by = 500))
{
  ind1 <- j-500+1
  nlsb1 <- nlsb[ind1:j]
  saveRDS(nlsb1, paste0("nlsb", ind1, "to", j, ".rds"))
  print(j)
}
nlsb1 <- nlsb[3501:3774]
saveRDS(nlsb1, paste0("nlsb", 3501, "to", 3774, ".rds"))

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
saveRDS(s0.true.xs, file = file.path(exp_dir, "s0_true_xs.rds")) 
saveRDS(s0.true.ra, file = file.path(exp_dir, "s0_true_ra.rds"))

# ------------------------------------------------------------------------------------------------------------
# Compute true A, WP

# Calculate true hydraulic parameters for cross sections
n.xs <- length(xWSEw)
A <- vector(length = n.xs)
WP <- vector(length = n.xs)
for (k in 1:n.xs)
{
  x <- cross_sections$x[[k]]
  b <- cross_sections$b[[k]]
  WSE <- max(b) # bankfull
  A[k] <- calc_A(x, b, WSE)
  WP[k] <- calc_WP(x, b, WSE)
}
saveRDS(A.true.xs, file = file.path(exp_dir, "A_true_xs.rds")) 
saveRDS(WP.true.xs, file = file.path(exp_dir, "WP_true_xs.rds")) 

# Calculate reach-average hydraulic parameters
n <- 2000 # number of segments in a reach = 10 km/5 m = 2000
A.r <- ra(A, n)
WP.r <- ra(WP, n)
saveRDS(A.r, file = file.path(exp_dir, "A_true_ra.rds")) 
saveRDS(WP.r, file = file.path(exp_dir, "WP_true_ra.rds")) 

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
saveRDS(w0.ra, file = "w0_ra.rds")
saveRDS(A0.true.ra, file = "A0_true_ra.rds")

# ------------------------------------------------------------------------------------------------
# Make plots of z0, A, WP along the river

k <- 12
plot(z0.true, main = paste("Minimum bed elevation (m) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(128,138))
lines(z0.l[,k], col = "red")
lines(z0.sb[,k], col = "orange")
lines(z0.sbm[,k], col = "purple")
lines(z0.nl[,k], col = "green")
lines(z0.nlsb[,k], col = "blue")
legend("bottomleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

plot(A.r, main = paste("Flow area (m^2) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(0,6000))
lines(A.l[,k], col = "red")
lines(A.sb[,k], col = "orange")
lines(A.sbm[,k], col = "purple")
lines(A.nl[,k], col = "green")
lines(A.nlsb[,k], col = "blue")
legend("bottomleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

plot(WP.r, main = paste("Wetted perimeter (m) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(450,850))
lines(WP.l[,k], col = "red")
lines(WP.sb[,k], col = "orange")
lines(WP.sbm[,k], col = "purple")
lines(WP.nl[,k], col = "green")
lines(WP.nlsb[,k], col = "blue")
legend("topleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# ------------------------------------------------------------------------------------------------
# Make plots of average z0, A, WP error at each exposure level

pred_z0 <- list(z0.l, z0.sb, z0.sbm, z0.nl, z0.nlsb)
pred_A <- list(A.l, A.sb, A.sbm, A.nl, A.nlsb)
pred_WP <- list(WP.l, WP.sb, WP.sbm, WP.nl, WP.nlsb)
plot_bias(expo, pred_z0, z0.true, na.rm = TRUE, 
          main = "z0 bias vs. exposure level, no meas. error", ylab = "Bias (m)", ylim = c(-2,2))
plot_bias(expo, pred_A, A.r, na.rm = TRUE, 
          main = "A bias vs. exposure level, no meas. error", ylab = "Bias (m)", ylim = c(-5000,6000))
plot_bias(expo, pred_WP, WP.r, na.rm = TRUE, 
          main = "WP bias vs. exposure level, no meas. error", ylab = "Bias (m)", ylim = c(-1,5))
legend("topright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# ------------------------------------------------------------------------------------------------
# Predict hydraulic parameters

# Initialize
z0.l <- array(dim = c(nr, n_exp_levels))
z0.sb <- array(dim = c(nr, n_exp_levels))
z0.sbm <- array(dim = c(nr, n_exp_levels))
z0.nl <- array(dim = c(nr, n_exp_levels))
z0.nlsb <- array(dim = c(nr, n_exp_levels))

A.l <- array(dim = c(nr, n_exp_levels))
A.sb <- array(dim = c(nr, n_exp_levels))
A.sbm <- array(dim = c(nr, n_exp_levels))
A.nl <- array(dim = c(nr, n_exp_levels))
A.nlsb <- array(dim = c(nr, n_exp_levels))

WP.l <- array(dim = c(nr, n_exp_levels))
WP.sb <- array(dim = c(nr, n_exp_levels))
WP.sbm <- array(dim = c(nr, n_exp_levels))
WP.nl <- array(dim = c(nr, n_exp_levels))
WP.nlsb <- array(dim = c(nr, n_exp_levels))

A0.l <- array(dim = c(nr, n_exp_levels))
A0.sb <- array(dim = c(nr, n_exp_levels))
A0.sbm <- array(dim = c(nr, n_exp_levels))
A0.nl <- array(dim = c(nr, n_exp_levels))
A0.nlsb <- array(dim = c(nr, n_exp_levels))

for (r in 1:nr) # takes about an hour to run the full loop
{
  for (k in 1:n_exp_levels)
  { 
    
    # if statements handle cases where the model is NULL/no model was fit
    if (class(lf[[r]][[k]]) == "lm")
    {
      z0.l[r,k] <- predict(lf[[r]][[k]], newdata = data.frame(w = 0))
      A.l[r,k] <- calc_model_A(lf[[r]][[k]], type = "linear")
      WP.l[r,k] <- calc_model_WP(lf[[r]][[k]], type = "linear")
      A0.l[r,k] <- calc_model_A0(lf[[r]][[k]], type = "linear")
    }
    if (class(sb[[r]][[k]][[1]]) == "lm")
    {
      z0.sb[r,k] <- predict(sb[[r]][[k]][[1]], newdata = data.frame(w = 0))
      A.sb[r,k] <- calc_model_A(sb[[r]][[k]], type = "sb")
      WP.sb[r,k] <- calc_model_WP(sb[[r]][[k]], type = "sb")
      A0.sb[r,k] <- calc_model_A0(sb[[r]][[k]], type = "sb")
    }
    if (any(class(sbm[[r]][[k]][[1]])=="lm"))
    {
      z0.sbm[r,k] <- predict(sbm[[r]][[k]][[1]], newdata = data.frame(w = 0))
      A.sbm[r,k] <- calc_model_A(sbm[[r]][[k]], type = "sbm")
      WP.sbm[r,k] <- calc_model_WP(sbm[[r]][[k]], type = "sbm")
      A0.sbm[r,k] <- calc_model_A0(sbm[[r]][[k]], type = "sbm")
    }
    if (class(nl[[r]][[k]]) == "nls")
    {
      z0.nl[r,k] <- predict(nl[[r]][[k]], newdata = data.frame(w = 0))
      A.nl[r,k] <- calc_model_A(nl[[r]][[k]], type = "nl", WSEw = rWSEw[[r]])
      WP.nl[r,k] <- calc_model_WP(nl[[r]][[k]], type = "nl", w = rWSEw[[r]]$w)
      A0.nl[r,k] <- calc_model_A0(nl[[r]][[k]], type = "nl", w0 = w0[r,k])
    }
    if (class(nlsb[[r]][[k]][[1]]) == "nls")
    {
      z0.nlsb[r,k] <- predict(nlsb[[r]][[k]][[1]], newdata = data.frame(w = 0))
      A.nlsb[r,k] <- calc_model_A(nlsb[[r]][[k]], type = "nlsb", WSEw = rWSEw[[r]]) # there may be a bug in the type = nlsb code here
      WP.nlsb[r,k] <- calc_model_WP(nlsb[[r]][[k]], type = "nlsb", w = rWSEw[[r]]$w)
      A0.nlsb[r,k] <- calc_model_A0(nlsb[[r]][[k]], type = "nlsb", w0 = w0[r,k])
    }
    
  }
  if (r%%10==0) {print(paste("progress:", r, "of", nr))}
}

# for (r in 1:nr) # just calculate area for nlsb method
# {
#   for (k in 1:n_exp_levels)
#   {
#     if (class(nlsb[[r]][[k]][[1]]) == "nls")
#     {
#       A.nlsb[r,k] <- calc_model_A(nlsb[[r]][[k]], type = "nlsb", WSEw = rWSEw[[r]])
#     }
#   }
#   print(r)
# }

save(z0.l, z0.sb, z0.sbm, z0.nl, z0.nlsb, file = file.path(exp_dir, "z0_pred.rda"))
save(A.l, A.sb, A.sbm, A.nl, A.nlsb, file = file.path(exp_dir, "A_pred.rda"))
save(WP.l, WP.sb, WP.sbm, WP.nl, WP.nlsb, file = file.path(exp_dir, "WP_pred.rda"))
save(A0.l, A0.sb, A0.sbm, A0.nl, A0.nlsb, file = file.path(exp_dir, "A0_pred.rda"))

# Compute slope via finite difference
s0.l <- diff(z0.l)
s0.sb <- diff(z0.sb)
s0.sbm <- diff(z0.sbm)
s0.nl <- diff(z0.nl)
s0.nlsb <- diff(z0.nlsb)
save(s0.l, s0.sb, s0.sbm, s0.nl, s0.nlsb, file = file.path(exp_dir, "s0_pred.rda"))

# ------------------------------------------------------------------------------------------------
# Calculate A0 prediction error for linear model with no measurement error

A0.var <- array(dim = c(nr, n_exp_levels))
for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    if (class(lf[[r]][[k]]) == "lm")
    {
      A0.var[r,k] <- calc_A0_variance(lf[[r]][[k]])
    }
  }
}
A0.sd <- sqrt(A0.var)
summary(A0.sd)
saveRDS(A0.sd, file.path(exp_dir, "A0_sd.rds"))

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

# ------------------------------------------------------------------------------------------------
# Plot z0, A0, S0, A, WP along the river at each exposure level

k <- 19 # choose an exposure level

# A0
plot(A0.true[,k], main = paste("A0 (m^2) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(0,500))
lines(A0.l[,k], col = "red")
lines(A0.sb[,k], col = "orange")
lines(A0.sbm[,k], col = "purple")
lines(A0.nl[,k], col = "green")
lines(A0.nlsb[,k], col = "blue")
legend("topleft", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# ------------------------------------------------------------------------------------------------
# Make plots of average z0, A, WP, s0 error at each exposure level

# z0
pred_z0 <- list(z0.l, z0.sb, z0.sbm, z0.nl, z0.nlsb)
z0.bias <- plot_bias(expo, pred_z0, z0.true, na.rm = TRUE,
          main = "z0 bias vs. exposure level, no meas. error", ylab = "Bias (m)", ylim = c(-25,3))
legend("bottomright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# s0
pred_s0 <- list(s0.l, s0.sb, s0.sbm, s0.nl, s0.nlsb)
# Convert units to cm/km
s0.true.cmkm <- lapply(s0.true, prod, 1e5)
pred.s0.cmkm <- lapply(pred_s0, prod, 1e5)
s0.bias <- plot_bias(expo, pred.s0.cmkm, s0.true.cmkm, na.rm = TRUE,
          main = "s0 bias vs. exposure level, no meas. error", ylab = "Bias (cm/km)")
legend("bottomright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# A
pred_A <- list(A.l, A.sb, A.sbm, A.nl, A.nlsb)
A.bias <- plot_bias(expo, pred_A, A.true, na.rm = TRUE,
          main = "A bias vs. exposure level, no meas. error", ylab = "Bias (m^2)")
legend("bottomright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# WP
pred_WP <- list(WP.l, WP.sb, WP.sbm, WP.nl, WP.nlsb)
WP.bias <- plot_bias(expo, pred_WP, WP.true, na.rm = TRUE,
          main = "WP bias vs. exposure level, no meas. error", ylab = "Bias (m)")
legend("bottomright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"),
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)

# A0
bias <- array(dim = c(n_exp_levels, 5)) # cannot use plot_bias because A0 changes with exposure level
na.rm = TRUE
for (k in 1:n_exp_levels)
{
  bias[k,1] <- mean(A0.l[,k] - A0.true[,k], na.rm = na.rm)
  bias[k,2] <- mean(A0.sb[,k] - A0.true[,k], na.rm = na.rm)
  bias[k,3] <- mean(A0.sbm[,k] - A0.true[,k], na.rm = na.rm)
  bias[k,4] <- mean(A0.nl[,k] - A0.true[,k], na.rm = na.rm)
  bias[k,5] <- mean(A0.nlsb[,k] - A0.true[,k], na.rm = na.rm)
}
bias <- as.data.frame(bias)
names(bias) <- c("l","sb","sbm","nl","nlsb")

# Plot bias vs. exposure level
plot(100*expo, bias$l, col = "red", type = "l", xlab = "Channel exposure (%)", 
     main = "A0 bias (m2), no meas. error", ylab = "bias (m)", xlim = c(40,100), ylim = c(0,1000))
lines(100*expo, bias$sb, col = "orange")
lines(100*expo, bias$sbm, col = "purple")
lines(100*expo, bias$nl, col = "green")
lines(100*expo, bias$nlsb, col = "blue")
abline(0,0)
legend("topright", legend = c("Zero", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1), ncol = 3)
saveRDS(bias, file = "A0_bias.rda")

# ------------------------------------------------------------------------------------------------





# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# --------------------------------------------- Scrap---------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Perform slope break method for all cross sections

M = 10
bias <- array(dim = c(nr, n_exp_levels))
z0.true <- vector(length = nr)
z0 <- array(dim = c(nr, n_exp_levels, M)) # reach, exposure, MC simulation iter
for (r in 1:nr)
{
  z0.true[r] <- rWSEw[[r]]$WSE[1] # rWSEw argument
  for (m in 1:n_exp_levels)
  {
    for (mm in 1:M)
    {
      WSEw_obs <- observe(rWSEw[[r]], exposure = expo[m])
      sbf <- fit_slopebreak(WSEw_obs, continuity = FALSE)
      if (!is.null(sbf))
      {
        z0[r,m,mm] <- coef(sbf[[1]])[1]
      }
    }
    bias[r,m] <- z0.bar[r,m] - z0.true[r]
  }
  print(paste0("Progress: ", 100*round(r/nr, 2), "%"))
}
z0.bar <- apply(z0, c(1,2), mean, na.rm = TRUE)
variance <- apply(z0, c(1,2), var, na.rm = TRUE)
# Started at 12:43 p.m. with 10 replicates for 3774 reach average cross sections
# Completed around 1:07 p.m. About 25 minutes for 10 replicates.
# 100 replicates will take around 4 hours?
save(z0, z0.true, z0.bar, bias, variance, file = file.path(fits_save_loc, sbsavename))


# Load models ------------------------------------------------------------------------------------

lf <- readRDS("/Users/jschap/Desktop/Cross_Sections/Outputs/nr3774_fitted_models_no_err/lf.rds")
sb <- readRDS("/Users/jschap/Desktop/Cross_Sections/Outputs/nr3774_fitted_models_no_err/sb.rds")
sbm <- readRDS("/Users/jschap/Desktop/Cross_Sections/Outputs/nr3774_fitted_models_no_err/sbm.rds")
nl1 <- readRDS("/Users/jschap/Desktop/Cross_Sections/Outputs/nr3774_fitted_models_no_err/nl_1_to_3500.rds")
nl2 <- readRDS("/Users/jschap/Desktop/Cross_Sections/Outputs/nr3774_fitted_models_no_err/nl3501to3774.rds")

vcov(nl2[[1]][[19]]) # variance of the estimated parameters

# Calculate prediction error

model <- lf[[1]][[19]]

# Using built-in R function
z0.1 <- predict(model, newdata = data.frame(w=0), se.fit = TRUE)
z0.1$se.fit # standard error of predicted mean

# Or manually-ish, assuming no measurement error:
cov.pars <- vcov(model)
sd_wse <- 0.1 # m
pred.var <- sd_wse^2 + cov.pars[1,1]

# Crashed at r=2724, k=14
# Error: at least one coef is NA: breakpoints at the boundary?


# ------------------------------------------------------------------------------------------------
# Model selection: Do split sample validation with each model

n_exp_levels <- 19
expo <- seq(0.05,0.95, by = 0.05)
sse <- array(dim = c(n_exp_levels, 5), data = 0)
for (k in 1:n_exp_levels)
{
  sum1 <- 0
  sum2 <- 0
  sum3 <- 0
  sum4 <- 0
  sum5 <- 0
  for (r in 1:nr)
  {
    sum1 <- sum1 + (A0.l[r,k] - A0.true[r,k])^2
    sum2 <- sum2 + (A0.sb[r,k] - A0.true[r,k])^2
    sum3 <- sum3 + (A0.sbm[r,k] - A0.true[r,k])^2
    sum4 <- sum4 + (A0.nl[r,k] - A0.true[r,k])^2
    sum5 <- sum5 + (A0.nlsb[r,k] - A0.true[r,k])^2
  }
  sse[k,1] <- sum1
  sse[k,2] <- sum2
  sse[k,3] <- sum3
  sse[k,4] <- sum4
  sse[k,5] <- sum5
}
sse <- as.data.frame(sse)
names(sse) <- c("l","sb","sbm","nl","nlsb")

plot(expo*100, sse[,1], 
     xlab = "channel exposure (%)", 
     ylab = "prediction error sse (m^4)", 
     type = "l", 
     col = "red",
     main = "prediction SSE for A0"
)
lines(expo*100, sse[,2], col = "orange")
lines(expo*100, sse[,3], col = "purple")
lines(expo*100, sse[,4], col = "green")
lines(expo*100, sse[,5], col = "blue")
legend("topright", legend = c("Linear","SB","SBM","NL","NLSB"), 
       col = c("red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)
saveRDS(sse, file = "A0_prediction_sse.rda")

# It's pretty clear that the NLSB method performs the best. However, there has been no accounting for model complexity.
# AIC and BIC can account for model complexity. However, calc_gof computes them wrt the fitted WSE-w values.
# I am interested in the goodness of fit for A0. Could AIC be calculated for this?

# ------------------------------------------------------------------------------------------------
# Plot A0 vs. predicted A0 at different exposure levels, and plot prediction error for linear models
# (This plot will look better if I do it in Matlab)

k <- 17
plot(A0.true[,k], 
     main = paste("A0 (m^2) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(0,500))
lines(A0.l[,k], col = "red")
lines(A0.l[,k] + 2*A0.sd[,k], col = "red", lty = 2)
lines(A0.l[,k] - 2*A0.sd[,k], col = "red", lty = 2)
legend("topleft", legend = "plus/minus 2 sd", col = "red", lty = 2)
