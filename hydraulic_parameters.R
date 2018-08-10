# Compute hydraulic parameters for each cross section

library(WSEw)

saveloc <- "/Users/jschap/Desktop/Cross_Sections/Data/Processed_Data"
load(file.path(saveloc, "processed_data_p21_sl_5m_hires.rda"))

# --------------------------------------------------------------------------------------------------

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

save(A, WP, file = file.path(saveloc, "xs_hydraul_params.rda"))

# Plot
par(mfrow = c(2,1))
plot(A, main = "Flow area (m^2)", type = "l")
plot(WP, main = "Wetted Perimeter (m)", type = "l")

# Calculate reach-average hydraulic parameters
n <- 2000 # number of segments in a reach = 10 km/5 m = 2000
A.r <- ra(A, n)
WP.r <- ra(WP, n)

# --------------------------------------------------------------------------------------------------

# Perform tests with idealized cross sections
load(file = "test_cross_sections.rda")

WSEw[[1]] # index 1 for rectangle, 2 for triangle, 3 for parabola

lf.test <- fit_linear(WSEw[[1]])
lf.test <- fit_linear(WSEw[[2]])
lf.test <- fit_linear(WSEw[[3]])

lf.test <- fit_nonlinear(WSEw[[1]])

calc_modelA()
calc_model_WP()


