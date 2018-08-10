# Executive file for WSEw project
# 
# Created 8/2/2018 JRS
# Revised 8/10/2018 JRS

# ------------------------------------------------------------------------------------------------

# Set up environment
library(raster)
library(strucchange)
library(segmented)
library(minpack.lm)
library(WSEw)

# Save names
fits_save_loc <- "/Users/jschap/Desktop/Cross_Sections/Outputs/Fitting_Results/"
sbsavename <- "r_nr3774_expo20_5m_0.05_sb.rda" # slope break
sbmsavename <- "r_nr3774_expo20_5m_0.05_sbm.rda" # multiple slope break
lsavename <- "r_nr3774_expo20_5m_0.05_l.rda" # linear
nlsavename <- "r_nr3774_expo20_5m_0.05_nl.rda" # nonlinear
nlsbsavename <- "r_nr3774_expo20_5m_0.05_nlsb.rda" # shape break

transects_name <- "/Users/jschap/Desktop/Cross_Sections/Data/Transects/p21_sl_5m_highres.rda"
processed_name <- "processed_data_p21_sl_5m_hires.rda" # savename for cross_sections, WSEw data

# ------------------------------------------------------------------------------------------------

# Load data

# Bathymetry
bathy.dir <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC"
bathy.name <- "bath_pool_21/bath_1999_p21"
bathyfile <- file.path(bathy.dir, bathy.name)
umesc <- raster(bathyfile)
levels(umesc)
depth_5 <- as_numeric_raster(umesc, att = "DEPTH_M") # native resolution depths (5 m)
writeRaster(depth_5, file = "Data/p21_depth.tif")
depth_5 <- raster("Data/p21_depth.tif")
depth_50 <- aggregate(depth_5, fact = 10) # resample depth to 50 m resolution
refWSE <- 470 # refWSE for pool 21 (ft)
refWSE <- refWSE/(39.37/12) #  convert to m
riv.dir <- "/Users/jschap/Desktop/Cross_Sections/Data/Centerlines"
load(file.path(riv.dir, "centerline21.rda"))
riv <- centerline_p21

# Processed cross section and WSE-w data
saveloc <- "/Users/jschap/Desktop/Cross_Sections/Data/Processed_Data"
load(file.path(saveloc, processed_name))

# ------------------------------------------------------------------------------------------------

# Process data prior to fitting

cross_sections <- auto_transects(section_length = 5, depth = depth_5, refWSE = refWSE, 
                                 savename = transects_name, makeplot = FALSE, riv = riv)
# (Takes about 4 hours at the highest possible resolution)

# Method 1: find corresponding flow width for WSE ranging from empty to bankfull conditions
xWSEw <- calc_WSEw(cross_sections, interval = 0.05, dx = 1) # number of data points depends on discretization
# Method 2: find corresponding WSE for flow width ranging from 100 m to bankfull width, in 50 m increments
# xWSEw <- calc_WSEw2(cross_sections, interval = 0.05, dx = 1) # anywhere from 7-21 data points, depending on river width

rWSEw <- reach_avg(xWSEw, l = 10000, res = 5)

save(cross_sections, xWSEw, rWSEw, file = file.path(saveloc, processed_name))

# ------------------------------------------------------------------------------------------------

# Fit linear model
lf <- fit_linear(WSEw_obs) # fit linear model

# Fit linear model, using Mersel method
lf.mersel <- fit_linear(WSEw_obs, mersel = TRUE, thres = 0.015) 

# Plot linear fit
plot(WSE~w, data = xWSEw[[1]], xlab = "width (m)", ylab = "WSE (m)",
     lty = 1, type = "l", lwd = 1, ylim = c(130,143))
points(WSE~w, data = WSEw_obs, col = "cyan")
lines(WSEw_obs$w, predict(lf), col = "blue", lwd = 1)
points(0, predict(lf, newdata = data.frame(w=0)), col = "blue", pch = 19, cex = 1)

# Calculate goodness of fit metrics
calc_gof(lf) # lapply it to all the linear fits

# --------------------------------------------------------------------------------------------------
# Fit slope break model

# 1. Linear with one breakpoint, continuous
sbf <- fit_slopebreak(WSEw_obs, continuity = TRUE)

# 2. Linear with one breakpoint, noncontinuous
sbf <- fit_slopebreak(WSEw_obs, continuity = FALSE)

# 3. Linear with several breakpoints, continous
sbf <- fit_slopebreak(WSEw_obs, continuity = TRUE, multiple_breaks = TRUE)

# 4. Linear with several breakpoints, noncontinuous (returns just the lowest fit)
sbf <- fit_slopebreak(WSEw_obs, multiple_breaks = TRUE, continuity = FALSE)

# 5. Mersel method: linear with one breakpoint, noncontinous,  screens out "non-optimal" cross sections
sbf <- fit_slopebreak(WSEw_obs, mersel = TRUE, thres = 0.015, window = 4)

sb.ind <- attributes(sbf)$sb.ind # get index of slope break

# Plot slope break fit
nn <- length(WSEw_obs$WSE)
plot(WSE~w, data = xWSEw[[1]])
points(WSEw_obs$w[sb.ind], WSEw_obs$WSE[sb.ind], pch = 19, col = "red")
lines(WSEw_obs$w[1:sb.ind], predict(sbf[[1]]))
lines(WSEw_obs$w[(sb.ind):nn], predict(sbf[[2]]))
points(0, predict(sbf[[1]], newdata = data.frame(w=0)), col = "blue", pch = 19, cex = 1)

# -------------------------------------------------------------------------------------------------
# Fit nonlinear model
fit <- fit_nonlinear(WSEw_obs)

# Plot nonlinear fit
plot(WSE~w, data = xWSEw[[1]], xlab = "width (m)", ylab = "WSE (m)",
     lty = 1, type = "l", lwd = 1, ylim = c(130,143))
points(WSE~w, data = WSEw_obs, col = "cyan")
lines(WSEw_obs$w, predict(fit), col = "blue", lwd = 1)
points(0, predict(fit, newdata = data.frame(w=0)), col = "blue", pch = 19, cex = 1)

# -------------------------------------------------------------------------------------------------
# Fit nonlinear slope break model

fits <- fit_nlsb(WSEw_obs)
sb.ind <- attributes(fits)$sb.ind # get index of slope break

# Plot nonlinear slope break fits
nn <- length(WSEw_obs$WSE)
plot(WSE~w, data = rWSEw[[1]])
points(WSEw_obs$w[sb.ind], WSEw_obs$WSE[sb.ind], pch = 19, col = "red")
lines(WSEw_obs$w[1:sb.ind], predict(fits[[1]]))
lines(WSEw_obs$w[(sb.ind):nn], predict(fits[[2]]))
points(0, predict(fits[[1]], newdata = data.frame(w=0)), col = "blue", pch = 19, cex = 1)

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Main experiments - model fitting at different exposure levels, known measurements
# Should use errors in variables method, but this can be refined later.
# Will do this with uncertain measurements later, preferably without MC simulation

expo <- seq(0.05, 0.95, by = 0.05) # exposure levels
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

begin.time <- Sys.time() # It does about 17 cross sections per minute.
for (r in 1:nr) # loop over reaches
{
  for (k in 1:n_exp_levels) # loop over exposure levels
  {
    WSEw_obs <- observe(WSEw = rWSEw[[r]], exposure = expo[k], sd_wse = 0, sd_w = 0)
    lf[[r]][[k]] <- fit_linear(WSEw_obs)
    sb[[r]][[k]] <- fit_slopebreak(WSEw_obs, multiple_breaks = FALSE, continuity = TRUE)
    sbm[[r]][[k]] <- fit_slopebreak(WSEw_obs, multiple_breaks = TRUE, continuity = TRUE)
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

# Crashed at r=2724, k=14
# Error: at least one coef is NA: breakpoints at the boundary?

save(lf, sb, sbm, nl, nlsb, file = file.path(saveloc, "nr10_fitted_models_no_err.rda"))

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

# ------------------------------------------------------------------------------------------------

# Evaluate performance of fits in terms of hydraulic parameters

load(file.path(saveloc, "r_hydraul_params_true.rda")) # load true hydraulic parameters
z0.true.xs <- unlist(lapply(cross_sections$b, min)) # get true bed elevations, though focus should be on hydraulic parameters, especially A0, WP0
z0.true <- ra(z0.true.xs, n = 2000)

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

# Compute z0, A, and WP predictions
for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  { 
    
    # if statements handle cases where the model is NULL/no model was fit
    if (class(lf[[r]][[k]]) == "lm")
    {
      z0.l[r,k] <- predict(lf[[r]][[k]], newdata = data.frame(w = 0))
      A.l <- calc_model_A(lf[[r]][[k]], type = "linear")
      WP.l <- calc_model_WP(lf[[r]][[k]], type = "linear")
    }
    if (class(sb[[r]][[k]][[1]]) == "lm")
    {
      z0.sb[r,k] <- predict(sb[[r]][[k]][[1]], newdata = data.frame(w = 0))
      A.sb <- calc_model_A(sb[[r]][[k]], type = "sb")
      WP.sb <- calc_model_WP(sb[[r]][[k]], type = "sb")
    }
    if (any(class(sbm[[r]][[k]][[1]])=="lm"))
    {
      z0.sbm[r,k] <- predict(sbm[[r]][[k]][[1]], newdata = data.frame(w = 0))
      A.sbm <- calc_model_A(sbm[[r]][[k]], type = "sbm")
      WP.sbm <- calc_model_WP(sbm[[r]][[k]], type = "sbm")
    }
    if (class(nl[[r]][[k]]) == "nls")
    {
      z0.nl[r,k] <- predict(nl[[r]][[k]], newdata = data.frame(w = 0))
      A.nl <- calc_model_A(nl[[r]][[k]], type = "nl", WSEw = rWSEw[[r]])
      WP.nl <- calc_model_WP(nl[[r]][[k]], type = "nl", w = rWSEw[[r]]$w)
    }
    if (class(nlsb[[r]][[k]][[1]]) == "nls")
    {
      z0.nlsb[r,k] <- predict(nlsb[[r]][[k]][[1]], newdata = data.frame(w = 0))
      A.nlsb <- calc_model_A(nlsb[[r]][[k]], type = "nlsb", WSEw = rWSEw[[r]])
      WP.nlsb <- calc_model_WP(nlsb[[r]][[k]], type = "nlsb", w = rWSEw[[r]]$w)
    }
    
  }
}

# Make plots
k <-5
plot(z0.true[1:10], main = paste("Minimum bed elevation (m) at", expo[k]*100,"% exposure"), 
     type = "l", ylim = c(127,149))
lines(z0.l[,k], col = "red")
lines(z0.sb[,k], col = "orange")
lines(z0.sbm[,k], col = "purple")
lines(z0.nl[,k], col = "green")
lines(z0.nlsb[,k], col = "blue")
legend("topright", legend = c("True", "Linear","SB","SBM","NL","NLSB"), 
       col = c("black", "red","orange", "purple","green","blue"), lwd = c(1,1,1,1,1,1))

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Scrap

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








