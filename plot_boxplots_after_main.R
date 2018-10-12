# This is just taken from the examples section of the parameter_predictions_boxplots documentation

# bplab <- c("1","1","1","1","2","2","2","2","3","3","3","3")
lumpedlab <- c("L","SB","NL","NLSB")
par(mfrow = c(2,2))
# par(opar)
k <- 4
# z0 error
parameter_predictions_boxplots(z0.l.error, z0.sb.error, z0.nl.error, z0.nlsb.error,
                               k = k,
                               lumped = TRUE,
                               main = paste("z0 error at", 100*expo[k], "percent channel exposure"),
                               ylab = "z0 error (m)",
                               legend = FALSE,
                               notch = TRUE,
                               names = lumpedlab,
                               ylim = c(-40,20))
# A0 error
k <- 16
parameter_predictions_boxplots(A0.l.error, A0.sb.error, A0.nl.error, A0.nlsb.error,
                               k = k,
                               lumped = TRUE,
                               main = paste("A0 error at", 100*expo[k], "percent channel exposure"),
                               ylab = "A0 error (sq. m)",
                               legend = FALSE,
                               notch = TRUE,
                               names = lumpedlab)

k <- 16
parameter_predictions_boxplots(z0.l, z0.sb, z0.nl, z0.nlsb, z0.true.ra, k = k, lumped = FALSE,
                               main = paste("z0 error at", 100*expo[k], "percent channel exposure"),
                               ylab = "average z0 error (m)", legend = FALSE, A0 = FALSE, notch = FALSE)
# for (r in 1:nr)
# {
#   abline(h = z0.true.ra[r])
# }

