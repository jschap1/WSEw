# This script is for determining what kind of river stages SWOT is likely to see
# Aug. 28, 2018 JRS

# Load gauge data
stage_name <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/miss_riv_at_lock_and_dam_21_2000-2018.txt"
gauge <- read.table(stage_name, header = TRUE, stringsAsFactors = FALSE)
gauge$Height <- as.numeric(gauge$Height)

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
h_bf <- 458.59+17
gauge <- read.table(stage_name, header = FALSE, stringsAsFactors = FALSE)
names(gauge) <- c("Month","Day","Year","Height")
gauge$Height <- as.numeric(gauge$Height)
h <- h[-which(h>h_bf)]
100*sum(h>h_bf)/length(h) # six percent of the stage values are flood flow


# Rescale and normalize before fitting, (also detrend)
h <- as.numeric(na.omit(gauge$Height))
z <- h/(max(h))

# data <- data.frame(t=1:length(z), z = z)
# m1 <- lm(z~t, data)
# z_dt <- residuals(m1) # detrended z

# z_rs <- (z_dt - min(z_dt))/(max(z_dt)-min(z_dt))
z_rs <- (z - min(z))/(max(z)-min(z))

z_rs[z_rs==0] <- 1e-3 # to avoid numerical errors

hist(z_rs, "fd")

fit1 <- fitdist(z_rs, distr = "burr", start = c(scale = 0.3, shape1 = 30, shape2 = 0.3))
fit1
dens1 <- dburr(seq(from = 0,to = 1,length.out = 500), 
      shape1 = as.numeric(fit1$estimate[2]),
      shape2 = as.numeric(fit1$estimate[3]),
      scale = as.numeric(fit1$estimate[1])
      )
#hist(z_rs, prob = TRUE, breaks = "fd", main = "Burr Fit", ylim = c(0,25))
#hist(z_rs, prob = TRUE, breaks = 100, main = "Burr Fit", ylim = c(0,30))
plot(density(z_rs))
lines(seq(from = 0,to = 1,length.out = 500), dens1, col="blue")

# Compare with gamma:
gammapar <- fitdistr(z_rs, "gamma")
g_shape <- as.numeric(gammapar$estimate[1])
g_rate <- as.numeric(gammapar$estimate[2])
g_dens <- dgamma(seq(from = 0,to = 1,length.out = 500), rate = g_rate, shape = g_shape)
lines(seq(from = 0,to = 1,length.out = 500), g_dens, col="darkred")

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
  a <- as.numeric(coef(m)[1]) # add trend back
  b <- as.numeric(coef(m)[2])
  z1 <- rep(a,n) + b*(1:n) + z_e
}

z.sim <- sim_z(scale, shape1, shape2, m1, min(z), max(z), n=100)
z.sim
summary(z.sim)

hist(z_rs1, "fd", prob=T)
hist(z_e, "fd", prob=T)

# --------------------------------------------------------------------------------------------
# What percent of bankfull depth will SWOT observe?
# Data-only approach, no fitting

stage_name <- "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/Stage/miss_at_quincy_1947-2018.txt"



which(h==min(h))

gauge <- read.table(stage_name, header = FALSE, stringsAsFactors = FALSE)
names(gauge) <- c("Month","Day","Year","Height")
gauge$Height <- as.numeric(gauge$Height)
h <- as.numeric(na.omit(gauge$Height))
summary(h)

d_bf <- 12.55 # m, from UMESC data near the Quincy gauge
h_datum <- 458.59 # from USACE gauge data website
flood_stage <- 17
h_bf <- h_datum + flood_stage
h_bf <- h_bf*0.3048
b.min <- h_bf - d_bf

# Calculate the fraction of bankfull depth time series
dbf_frac <- (h*0.3048 - b.min)/d_bf

plot(dbf_frac, main = "Fraction of bankfull depth", type = "l")

# OK, now pretend you are SWOT:
ndays <- 10
T <- length(dbf_frac)
swot_sample <- vector(length = ndays, "list")
start_day <- 1
while (start_day < ndays)
{
  swot_sample[[start_day]] <- dbf_frac[seq(start_day, T, by = ndays)]
  start_day <- start_day + 1
}

hist(dbf_frac, "fd")
k <- 3
hist(swot_sample[[k]], "fd")


