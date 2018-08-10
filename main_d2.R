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
fits_save_loc <- "/Users/jschap/Desktop/Cross_Sections/Data/Fitting_Results/w3"
sbsavename <- "r_nr3774_expo20_5m_0.05_sb.rda" # slope break
lsavename <- "r_nr3774_expo20_5m_0.05_l.rda" # linear
nlsavename <- "r_nr3774_expo20_5m_0.05_nl.rda" # nonlinear
nlsbsavename <- "r_nr3774_expo20_5m_0.05_nlsb.rda" # shape break

savename <- "/Users/jschap/Desktop/Cross_Sections/Data/Transects/p21_sl_5m_highres.rda"

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
load(file.path(saveloc, "processed_data_p21_sl_5m_hires.rda"))

# ------------------------------------------------------------------------------------------------

# Process data prior to fitting

cross_sections <- auto_transects(section_length = 5, depth = depth_5, refWSE = refWSE, 
                                 savename = savename, makeplot = FALSE, riv = riv)
# (Takes about 4 hours at the highest possible resolution)

# Method 1: find corresponding flow width for WSE ranging from empty to bankfull conditions
xWSEw <- calc_WSEw(cross_sections, interval = 0.05, dx = 1) # number of data points depends on discretization
# Method 2: find corresponding WSE for flow width ranging from 100 m to bankfull width, in 50 m increments
# xWSEw <- calc_WSEw2(cross_sections, interval = 0.05, dx = 1) # anywhere from 7-21 data points, depending on river width

rWSEw <- reach_avg(xWSEw, l = 10000, res = 5)

save(cross_sections, xWSEw, rWSEw, file = file.path(saveloc, "processed_data_p21_sl_5m_hires.rda"))

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


