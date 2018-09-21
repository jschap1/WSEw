setwd("/Users/jschap/Desktop/Cross_Sections")
opar <- par()
load("Data/Processed/processed_xs_data.rda")

r <- 3000
par(mfrow = c(1,2))


# original
plot(cross_sections$x[r][[1]], cross_sections$b[r][[1]], type="l",
     xlab = "x (m)", ylab = "b (m)", main = "Example cross section")
plot(WSE~w, xWSEw[[r]], xlab = "width (m)", ylab = "height (m)", main = "Height-width samples")


# reach averaged (smoothed)
plot(WSE~w, rWSEw[[r]])