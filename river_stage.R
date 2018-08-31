# This script is for determining what kind of river stages SWOT is likely to see
# Aug. 28, 2018 JRS

# Load gauge data
stage_name <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/miss_riv_at_lock_and_dam_21_2000-2018.txt"
gauge <- read.table(stage_name, header = TRUE, stringsAsFactors = FALSE)
gauge$Height <- as.numeric(gauge$Height)

# Add gauge location to an existing study area plot
gauge_coords <- data.frame(lon = -91.42904167, lat = 39.90380833)
coordinates(gauge_coords) <- ~lon+lat
crs(gauge_coords) <- "+init=epsg:4326"
gauge_coords.utm <- spTransform(gauge_coords, crs(riv))
points(gauge_coords.utm, pch = 19, col = "blue")

# Sample gauge data as SWOT would (every ~10 days)
# Use all possible start date/combos

N <- dim(gauge)[1]
swot_inds <- vector(length = length(ind1), "list")
ind1 <- 1:10
for (k in 1:length(ind1))
{
  swot_inds[[k]] <- seq(ind1[k], N, by = 10)  
}

# Plot histogram of all the data
hist(gauge$Height, breaks = "FD")
# hist(gauge$Height, breaks = "Sturges")

# Plot histogram of each of the possible SWOT samplings
# Going to skip this for now and treat the whole set of heights

# Plot probability density of the data
plot(density(gauge$Height, na.rm = TRUE), 
     main = "Empirical PDF for gauge data", 
     xlab = "Stage (ft)"
     )

# Repeat for logs
plot(density(log(gauge$Height), na.rm = TRUE), 
     main = "Empirical PDF for log gauge data", 
     xlab = "Stage (ft)"
     ) # does not make it normal, at all, so this rules out lognormal

# Fit gamma density
library(MASS)
x <- gauge$Height
x <- na.omit(x)
fit.params <- fitdistr(x, "gamma", lower = c(0,0))
summary(fit.params)

y.gamma <- dgamma(x = seq(min(x),max(x),length.out = 500), 
       shape = fit.params$estimate[1],
       rate = fit.params$estimate[2]
       )
# Plot the gamma density
lines(x = seq(min(x),max(x),length.out = 500), y, col = "red")

# Fit Weibull density
fit.params <- fitdistr(x, "weibull", lower = c(0,0))
y.weibull <- dweibull(x = seq(min(x),max(x),length.out = 500), 
            shape = fit.params$estimate[1],
            scale = fit.params$estimate[2]
)
lines(x = seq(min(x),max(x),length.out = 500), y.weibull, col = "green")

# Fit lognormal density
fit.params <- fitdistr(x, "lognormal", lower = c(0,0))
y.lnorm <- dlnorm(x = seq(min(x),max(x),length.out = 500), 
                      meanlog = fit.params$estimate[1],
                      sdlog = fit.params$estimate[2]
)
lines(x = seq(min(x),max(x),length.out = 500), y.lnorm, col = "yellow")

# Fit exponential density (account for shift)
x.shift <- x-min(x)
fit.params <- fitdistr(x.shift, "exponential")
y.exp <- dexp(x = seq(min(x.shift),max(x.shift),length.out = 500), 
                  rate = fit.params$estimate[1]
)
lines(x = seq(min(min(x)+x.shift),max(min(x) + x.shift),length.out = 500), y.exp, col = "blue")

legend("topright", 
       legend = c("empirical","gamma", "weibull", "lognormal", "exponential"), 
       fill = c("black","red","green", "yellow", "blue")
       )

# ----------------------------------------------------------------------------------------------------
# Perform final fit

WSE <- na.omit(gauge$Height)
WSEmax <- max(WSE)
x <- WSE/WSEmax
x.shift <- x-min(x)
fit.params <- fitdistr(x.shift, "exponential")

plot(density(x, na.rm = TRUE))

y.exp <- dexp(x = seq(min(x.shift),max(x.shift),length.out = 500), 
              rate = fit.params$estimate[1]
)

lines(x = seq(min(min(x)+x.shift),max(min(x) + x.shift),length.out = 500), y.exp, col = "blue")

# Use the PDF of the gauge heights to randomly sample the possible WSE values as if making SWOT measurements

# ----------------------------------------------------------------------------------------------------
# Perform final final fit (turns out exponential is problematic)

fit.params <- fitdistr(x, "gamma")
plot(density(x, na.rm = TRUE))
y.gamma <- dgamma(x = seq(min(x)-0.1,max(x),length.out = 500),
                  shape =  fit.params$estimate[1], 
                  rate =  fit.params$estimate[2]
                  )
lines(seq(min(x)-0.1,max(x),length.out = 500), y.gamma, col = "red")

# What about the Gumbel distribution? ("Generalized extreme value distribution type I)
# It is used to model annual peak flows. Look at some papers on the subject of 
# "Probability distributions for daily river stage."

# ----------------------------------------------------------------------------------------------------
# Burr distribution workflow

library(fitdistrplus)
library(actuar)
setwd("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB")

# Load gauge data
stage_name <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/Stage/miss_at_quincy_1947-2018.txt"
# stage_name <- "/Data/Stage/miss_at_quincy_1947-2018.txt"
# stage_name <- "/Data/miss_riv_at_lock_and_dam_21_2000-2018.txt"
gauge <- read.table(stage_name, header = FALSE, stringsAsFactors = FALSE)
names(gauge) <- c("Month","Day","Year","Height")
gauge$Height <- as.numeric(gauge$Height)

# Rescale and normalize before fitting, also detrend
h <- as.numeric(na.omit(gauge$Height))
h <- na.omit(gauge$Height)
hbf <- max(h)
hbf <- 478.5
z <- h/hbf

z_rs <- (z - min(z))/(max(z)-min(z))
hist(z_rs, "fd")

data <- data.frame(t=1:length(z), z = z)
m1 <- lm(z~t, data)
z_dt <- residuals(m1) # detrended z

zmin <- min(z)
zmax <- max(z)
z_rs <- (z_dt - min(z_dt))/(max(z_dt)-min(z_dt))
# z_rs <- (z - min(z))/(max(z)-min(z))

z_rs[z_rs==0] <- 1e-3 # to avoid numerical errors

hist(z_rs, "fd")

fit1 <- fitdist(as.numeric(z_rs), distr = "burr", start = c(scale = 0.3, shape1 = 30, shape2 = 0.3))
fit1
dens1 <- dburr(seq(from = 0,to = 1,length.out = 500), 
      shape1 = as.numeric(fit1$estimate[2]),
      shape2 = as.numeric(fit1$estimate[3]),
      scale = as.numeric(fit1$estimate[1])
      )
hist(z_rs, prob = TRUE, breaks = "fd", main = "Burr Fit", ylim = c(0,25))
hist(z_rs, prob = TRUE, breaks = 100, main = "Burr Fit", ylim = c(0,30))

plot(density(z_rs), main = "Burr fit")
lines(seq(from = 0,to = 1,length.out = 500), dens1, col="blue")
legend("topright", legend = c("Observed","Fitted"), lty = c(1,1), col = c("black","blue"))

#' Simulate SWOT WSE
#'
#' Simulate z values for SWOT sampling using the fitted Burr distribution
#' @export
#' @param scale
#' @param shape1
#' @param shape2
#' @param m linear model used to detrend the data
#' @param n number of replicates
#' @example z.sim <- sim_z(scale, shape1, shape2, m1, min(z), max(z), n=100)

shape1 = as.numeric(fit1$estimate[2])
shape2 = as.numeric(fit1$estimate[3])
scale = as.numeric(fit1$estimate[1])

sim_z <- function(scale, shape1, shape2, m, zmin, zmax, n)
{
  z_rs1 <- rburr(n, shape1 = shape1, shape2 = shape2, scale = scale)
  z_e <- zmin + z_rs1*(zmax - zmin) # scale, called z_e because it is the residual from the regression of z on time
#  a <- as.numeric(coef(m)[1]) # add trend back
 # b <- as.numeric(coef(m)[2])
  #z1 <- rep(a,n) + b*(1:n) + z_e
}

z.sim <- sim_z(scale, shape1, shape2, 1, zmin, zmax, n=20000)
z.sim
summary(z.sim)
summary(z)

par(mfrow = c(2,1))
hist(z.sim, "fd", prob=T, xlim =c(0.97, 1.02))
hist(z, "fd", prob=T, xlim =c(0.97, 1.02))

# ------ Figure for the Paper
# Simulate some SWOT measurements using the fitted Burr distribution

z.sim <- sim_z(scale, shape1, shape2, 1, zmin, zmax, n=100)
hbf.xs <- max(xWSEw[[3000]]$WSE)
sim.obs <- z.sim*hbf.xs # simulated SWOT observations
xWSEw3 <- calc_WSEw3(cross_sections, dist = "burr", n.obs = floor(1*365/10)) # 3 yrs of obs

plot(WSE~w, xWSEw[[3000]], type = "l", main = "Example h-w relationship", 
     xlab = "width (m)", ylab = "h (m)")

points(w.obs, sim.obs)
points(WSE~w, xWSEw3[[3000]], col = "red", pch=19)
legend("topleft", 
       legend = c("h-w relationship","sampled h-w values"), 
       col = c("black", "red"), 
       lty = c(1, NA), pch = c(NA, 19)
       )

# ------------------------------------------------------------------------------------------------------------------------------
# Get the 2-year flood stage

# Load stage data
stage_name <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/Stage/miss_at_quincy_1947-2018.txt"
gauge <- read.table(stage_name, header = FALSE, stringsAsFactors = FALSE)
names(gauge) <- c("Month","Day","Year","Height")
gauge$Height <- as.numeric(na.omit(gauge$Height))
hist(gauge$Height)
head(gauge)
tail(gauge)
timevector <- seq(as.Date("1947-01-01"), as.Date("2018-08-29"), by = "1 day")

# remove the empty times

# check length is the same
length(timevector)
dim(gauge)[1]

save(timevector, gauge, file = "Data/quincy_gauge_data.rda")

# Add water year column to df

wtr_yr <- function(dates, start_month = 10)
{ # credit to Nick Rong at https://blogs.ubc.ca/nickrong/2015/09/22/simple-function-in-r-to-calculate-water-year
  # Year offset
  offset = ifelse(as.integer(format(dates, "%m")) < start_month, 0, 1)
  # Water year
  adj.year = as.integer(format(dates, "%Y")) + offset
  # Return the water year
  return(adj.year)
}

gauge$wtr_yr <- wtr_yr(timevector)

# Get peak annual flows
wy <- unique(gauge$wtr_yr)
n.wy <- length(wy)
h.peak <- vector(length = n.wy)
for (k in 1:n.wy)
{
  h.peak[k] <- max(gauge$Height[gauge$wtr_yr == wy[k]], na.rm = TRUE)
}

Fn <- ecdf(h.peak)
Fn(478.5) # this is approximately the 2-year peak flood stage

plot(ecdf(h.peak))


