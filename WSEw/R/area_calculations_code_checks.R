# Load test cross sections

library(WSEw)
source("./Codes/test_case.R")

# Triangular cross section

r <- 2
plot(cross_sections$x[[r]], cross_sections$b[[r]], type = "l")

trueA <- 0.5*25*100
trueWP <- 2*sqrt(50^2+25^2)

xWSEw <- calc_WSEw(cross_sections, interval = 0.05, dx = 1)

# Calculate hydraulic parameters using different methods

A_from_WSEw <- calc_A_from_WSEw(xWSEw[[r]])

A_from_calc_A <- calc_A(cross_sections$x[[r]], cross_sections$b[[r]], WSE = 25)

m <- fit_linear(xWSEw[[r]])
A_from_calc_model_A0 <- calc_model_A0(model = m, type = "linear")

WP <- calc_WP(cross_sections$x[[r]], cross_sections$b[[r]], WSE = 25)
WP_model <- calc_model_WP(model = m, type = "linear")

# Check that they are correct

# See what happens when you add error       
