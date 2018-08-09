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
