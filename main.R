# Executive file for WSEw project
# 
# Created 8/2/2018 JRS

# --------------------------------------------------------------------------------------------------
# Set up environment

rm(list=ls())

# Ideally, these packages would be imported by the WSEw package.
# Better yet, it wouldn't require this many libraries; it is a lot.
# library(rgdal) 
# library(maptools)
library(raster)
library(strucchange)
library(segmented)
library(minpack.lm)
library(WSEw)

setwd("/Users/jschap/Desktop/Cross_Sections")

# Save names
fits_save_loc <- "/Users/jschap/Desktop/Cross_Sections/Data/Fitting_Results/w3"
sbsavename <- "r_nr3774_expo20_5m_0.05_sb.rda" # slope break
lsavename <- "r_nr3774_expo20_5m_0.05_l.rda" # linear
nlsavename <- "r_nr3774_expo20_5m_0.05_nl.rda" # nonlinear
nlsbsavename <- "r_nr3774_expo20_5m_0.05_nlsb.rda" # shape break

# --------------------------------------------------------------------------------------------------
# Load data

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

savename <- "/Users/jschap/Desktop/Cross_Sections/Data/Transects/p21_sl_5m_highres.rda"

riv.dir <- "/Users/jschap/Desktop/Cross_Sections/Data/Centerlines"
load(file.path(riv.dir, "centerline21.rda"))
riv <- centerline_p21

# --------------------------------------------------------------------------------------------------
# Process data prior to fitting

cross_sections <- auto_transects(section_length = 5, depth = depth_5, refWSE = refWSE, 
               savename = savename, makeplot = FALSE, riv = riv)
# (Takes about 4 hours at the highest possible resolution)

# Method 1: find corresponding flow width for WSE ranging from empty to bankfull conditions
xWSEw <- calc_WSEw(cross_sections, interval = 0.05, dx = 1) # number of data points depends on discretization
# Method 2: find corresponding WSE for flow width ranging from 100 m to bankfull width, in 50 m increments
# xWSEw <- calc_WSEw2(cross_sections, interval = 0.05, dx = 1) # anywhere from 7-21 data points, depending on river width

rWSEw <- reach_avg(xWSEw, l = 10000, res = 5)

saveloc <- "/Users/jschap/Desktop/Cross_Sections/Data/Processed_Data"
load(file.path(saveloc, "processed_data_p21_sl_5m_hires.rda"))
save(cross_sections, xWSEw, rWSEw, file = file.path(saveloc, "processed_data_p21_sl_5m_hires.rda"))

# Make observations
WSEw_obs <- observe(xWSEw[[1]], exposure = 1, sd_wse = 0, sd_w = 0)

# --------------------------------------------------------------------------------------------------
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

# -------------------------------------------------------------------------------------------------
# Test modeling across all different exposure levels, but with no measurement error
# Can do for reach-average or individual cross sections.

lsavename <- "r_nr3774_expo20_5m_0.05_l.rda" # linear

nr <- length(rWSEw) # number of reaches
expo <- seq(0.05, 0.95, by = 0.05)
n_exp_levels <- length(expo)
bias <- array(dim = c(nr, n_exp_levels))
z0.true <- vector(length = nr)
z0<- array(dim = c(nr, n_exp_levels))
begin.time <- Sys.time()
for (r in 1:nr)
{
  for (m in 1:n_exp_levels)
  {
    WSEw_obs <- observe(rWSEw[[1]], exposure = 1, sd_wse = 0, sd_w = 0) #!
    lf <- fit_linear(WSEw_obs)
    z0[r,m] <- coef(lf[1])[1]
    z0.true[r] <- rWSEw[[r]]$WSE[1] #!
    bias[r,m] <- z0[r,m] - z0.true[r]
  }
  print(paste0("Progress: ", 100*round(r/nr, 2), "%"))
}
duration <- Sys.time() - begin.time
print(paste("Took", round(duration, 2), "minutes")) # 2 minutes to run for all reach-avg cross sections
save(z0, z0.true, bias, file = file.path(fits_save_loc, lsavename))


# -------------------------------------------------------------------------------------------------
# Slope break method with measurement error at different exposure levels

# Perform slope break method for all cross sections
n_exp_levels <- length(expo)
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

# -------------------------------------------------------------------------------------------------
# Commands for updating the WSEw R package

library(devtools)
library(roxygen2)
setwd("Codes/pkg")
rm(list=ls())
getwd()
setwd("pkg/WSEw")
setwd("..")
document()
install("WSEw")
library(WSEw)
?resample_polyline
?auto_transects
resample_polyline
check()
