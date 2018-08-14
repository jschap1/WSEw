# WSEw project - code testing
# 
# Created 8/7/2018 JRS

# --------------------------------------------------------------------------------------------------
# Set up environment and load data

rm(list=ls())
#library(raster)
library(WSEw)
setwd("/Users/jschap/Desktop/Cross_Sections")

bathy.dir <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC"
bathy.name <- "bath_pool_21/bath_1999_p21"
bathyfile <- file.path(bathy.dir, bathy.name)

umesc <- raster(bathyfile)
levels(umesc)
depth_5 <- raster("Data/p21_depth.tif")

refWSE <- 470 # refWSE for pool 21 (ft)
refWSE <- refWSE/(39.37/12) #  convert to m

riv.dir <- "/Users/jschap/Desktop/Cross_Sections/Data/Centerlines"
load(file.path(riv.dir, "centerline21.rda"))
riv <- centerline_p21

saveloc <- "/Users/jschap/Desktop/Cross_Sections/Data/Processed_Data"
load(file.path(saveloc, "processed_data_p21_sl_5m_hires.rda"))
save(cross.sections, file = file.path(saveloc, "cross_sections_format2.rda"))

# -------------------------------------------------------------------------------------------------
# Simple cross section

x <- seq(0,600,100)
b <- c(144,136,135,136,135,136,144)

cross.section <- list(x = x, b = b, d = b) 

# WSEw <- calc_WSEw(cross.section, interval = 0.05, dx = 1)
# Calculate it using calc_WSEw and get_width. Must be run manually.

xn <- 1000
r <- 1
cross.section <- list(x = cross_sections$x[[xn]], b = cross_sections$b[[xn]], d = cross_sections$d[[xn]])
WSEw <- xWSEw[[xn]]

par(mfrow = c(1,2))
plot(b~x, cross.section, type = "l")
plot(WSE~w, WSEw)

lf1 <- fit_linear(WSEw)
sb1 <- fit_slopebreak(WSEw, multiple_breaks = FALSE, continuity = TRUE)
sbm1 <- fit_slopebreak(WSEw, multiple_breaks = TRUE, continuity = TRUE)
nl1 <- fit_nonlinear(WSEw)
nlsb1 <- fit_nlsb(WSEw)

A.true <- calc_A(cross.section$x, cross.section$b, WSE = max(WSEw$WSE))
A.true <- calc_A_from_WSEw(WSEw)
A.l <- calc_model_A(lf1, type = "linear")
A.sb <- calc_model_A(sb1, type = "sb")
A.sbm <- calc_model_A(sbm1, type = "sbm")
A.nl <- calc_model_A(nl1, type = "nl", WSEw)
A.nlsb <- calc_model_A(nlsb1, type = "nlsb", WSEw)

z0.l <- predict(lf1, newdata = data.frame(w=0))
z0.sb <- predict(sb1[[1]], newdata = data.frame(w=0))
z0.sbm <- predict(sbm1[[1]], newdata = data.frame(w=0))
z0.nl <- predict(nl1, newdata = data.frame(w=0))
z0.nlsb <- predict(nlsb1[[1]], newdata = data.frame(w=0))

# Plot modeled cross sections over true cross section
par(mfrow = c(1,1))
plot(b~x, cross.section, type = "l", ylim = c(130, 145), main = "Simple cross section test", lwd = 2)
plot_model(lf1, type = "linear", col = "red")
plot_model(sb1, type = "sb", col = "purple")
plot_model(nl1, type = "nl", WSEw, col = "blue")
plot_model(nlsb1, type = "nlsb", WSEw, col = "green")
legend("top", ncol = 2, legend = c("true","linear","sb","nl", "nlsb"), fill = c("black","red","purple","blue", "green"))

# -------------------------------------------------------------------------------------------------
# calc_WSEw / get_width

#cross.section <- list(x = cross_sections$x[[1]], b = cross_sections$b[[1]], d = cross_sections$d[[1]]) 
# may want to update to a more useful format, instead of a list of length 3, better to be a list of length nseg,
# except this would interfere with other functions that are set up for cross_sections to be in a particular format

# Using the updated format (format2) for the cross sections:
#x <- cross.sections[[1]]$x
#b <- cross.sections[[1]]$b
#d <- cross.sections[[1]]$d

#WSEw <- calc_WSEw(cross.section, interval = 0.05, dx = 1)




