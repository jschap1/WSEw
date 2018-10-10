# Executive file for WSEw project, but for the Pepsi Challege data

# ------------------------------------------------------------------------------------------------
# Set up environment

library(foreach)
library(tools)
library(doMC)
library(raster)
library(strucchange)
library(minpack.lm)
library(WSEw)

# library(devtools)
# library(roxygen2)
# update_WSEw()

setwd("/Users/jschap/Documents/Research/SWOTBATH")

opar <- par()

# Make a directory to store results
exp_desc <- "garonne_10km_spacing_5m_sampling_even_err_type_MC_replicates_500"
exp_dir <- file.path("./Outputs", exp_desc) # directory for this experiment's outputs

M <- 500 # number of MC replicates

expo <- seq(0.05, 0.95, length.out = 19) # exposure levels
n_exp_levels <- length(expo)

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

rivername <- "GaronneDownstream_rWSEw.rds"
rWSEw <- readRDS(file.path("./Data/Pepsi_Challenge_I/rWSEw_Pepsi1", rivername))

nr <- length(rWSEw)

truth_dir <- file.path("./true_parameters", strsplit(file_path_sans_ext(rivername), "_rWSEw")[[1]])
if(!dir.exists(truth_dir))
{
  dir.create()
}

par(mfrow = c(3,3))
for (r in 1:nr)
{
  plot(WSE~w, rWSEw[[r]], xlab = "width (m)", ylab = "height (m)", 
       main = "Height-width relationship", type = "o")
}

# ------------------------------------------------------------------------------------------------
# Calculate true parameters

# z0
z0.true.ra <- vector(length = nr)
for (r in 1:nr)
{
  z0.true.ra[r] <- min(rWSEw[[r]]$WSE)
}
saveRDS(z0.true.ra, file = file.path(truth_dir, "z0_true_ra.rds"))

# s0
s0.true.ra <- diff(z0.true.ra)
saveRDS(s0.true.ra, file = file.path(truth_dir, "s0_true_ra.rds"))

# par(mfrow=c(2,1))
# plot(s0.true.ra, type = "l")
# abline(0,0)

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
# Takes about 15 minutes for 8 reach average cross sections (Garonne Downstream)

# Run the parallel computations (see computation_time_calculator.xls)

# the models are much larger and slower to fit bc there are many more data points (not nn = 21, anymore) 
# and also maybe because the rows are named
# to make them faster to run and to make them more comparable with the UMRB case, I will thin out the data,
# using even sampling ********

tWSEw <- thin_evenly(rWSEw[[1]], nn = 21) # use a function similar to the one I used to rescale the cross sections prior to averaging

resample_xs(cross_sections, n = 300) # Adapt resample_xs to do this.


ncores <- detectCores()
registerDoMC(cores = max(nr, ncores - 1))

begin.time <- Sys.time()
WSEw_val <- foreach(r = 1:nr, .combine = c) %dopar% {observe_par(r, rWSEw)} # returns a dummy value; the point is to save .rds files
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

# Takes about 5 min. for 8 reaches (Garonne Downstream)

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
# This is physically inconsistent, so can neglect these predictions and analyze the remaining predicted values.

# pred_lf <- clean_by_A0(pred_lf, h1) # remove A0 and z0 predictions where A0<0
# pred_sb <- clean_by_A0(pred_sb, h1)
# pred_nl <- clean_by_A0(pred_nl, h1)
# pred_nlsb <- clean_by_A0(pred_nlsb, h1)

# Run predictions again (if necessary) for any reaches where errors occurred

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

# ------------------------------------------------------------------------------------------------
# Load data from previous runs (optional)

# Cross section and WSE-w data
load(file.path(truth_dir, "processed_xs_data_10km.rda"))
# loads {cross_sections, cross_sections_avg, xWSEw, rWSEw}

# Fitted models
lf <- readRDS(file.path(exp_dir, "lf/lf_r_1.rds"))
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
