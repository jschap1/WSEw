# Fit an ARIMA model to the river stage data
# Code based on STAT221 project from Winter 2018

setwd("/Users/jschap/Documents/Research/SWOTBATH/")

stage <- read.table("./Data/USACE/Stage/miss_riv_at_lock_and_dam_21_2000-2018.txt", 
                    header=TRUE, stringsAsFactors = FALSE)

stage <- read.table("Stage/miss_at_quincy_1947-2018.txt", 
                    stringsAsFactors = FALSE, header = TRUE)

stage$Height <- as.numeric(stage$Height)
n <- dim(stage)[1]
head(stage)
tail(stage)
# t <- seq(as.Date("2000-01-01"), as.Date("2018-08-28"), by = "1 day")
t <- seq(as.Date("1947-01-01"), as.Date("2018-08-29"), by = "1 day")

# Find missing dates
missing <- which(diff(stage$Day)>1)
k <- 1
stage[(missing[k]-3):(missing[k]+3),]

stage[3703:3709,]

# Missing dates (Lock and dam 21: 
# 11/15/2005
# 7/21/2008
# 7/22/2008
# 4/4/2009
# 2/27/2010
# 2/28/2010
# 6/29/2014
# 1/18/2015
# 11/27/2015

# Missing dates (Quincy): 
# 10/21/2011
# 10/22/2011
# 10/23/2011

# Eliminate the times corresponding to missing dates
elim.ind <- which(t == "2005-11-15" | t == "2008-07-21" | t == "2008-07-22" | t == "2009-04-04"
      | t == "2014-06-29" | t == "2015-01-18" | t == "2015-11-27" | t == "2010-02-27" | t == "2010-02-28")
t <- t[-elim.ind]

elim.ind <- which(t == "2010-10-21" | t == "2010-10-22" | t == "2010-10-23")
t <- t[-elim.ind]

# Save the stage and time data
h <- stage$Height
h <- h*0.3048 # convert to meters
na.ind <- which(is.na(h)) # omit missing values
t <- t[-na.ind]
h <- na.omit(h)
save(t, h, file = "stage_height_quincy.rda")

# -----------------------------------------------------------
# Exploratory data analysis

load("stage_height.rda")

par(mfrow = c(2,2))

plot(t, h, main = "Stage (m) vs. time", type = "l")

hist(h, ann = FALSE, freq = FALSE)
x = seq(min(h),max(h),length=1000)
lines(x, dnorm(x, mean = mean(h), sd = sd(h)), type="l", col = "blue")
legend(locator(1), c("Stage","Normal Distribution"), lty = c(1,1), col = c("black", "blue"))
title("PDF of daily river stage")

plot(t, log(h), main = "log Stage vs. time", type = "l")

hist(log(h), ann = FALSE, freq = FALSE)
x = seq(min(log(h)),max(log(h)),length=1000)
lines(x, dnorm(x, mean = mean(log(h)), sd = sd(log(h))), type="l", col = "blue")
title("PDF of monthly log stage")

# ------------------------------------------------------------------------------------
# Plot ACFs and PACFs

acf.la <- acf(log(h), lag.max = 100)
pacf.la <- pacf(log(h), lag.max = 100)

par(mfrow = c(2,1))
plot(acf.la$lag[-1], acf.la$acf[-1], type = "h", main = "ACF, log h",
     xlab = "Lag", ylab = "ACF")
lines(x = c(0,acf.la$lag[101]), y = c(0,0), lty = 1)
plot(pacf.la$lag[-1], pacf.la$acf[-1], type = "h", main = "PACF, log h",
     xlab = "Lag", ylab = "PACF")
lines(x = c(0,length(pacf.la$lag)), y = c(0,0), lty = 1)

# ------------------------------------------------------------------------------------
# Need to de-trend/de-seasonalize the data

n <- length(t)

df <- as.data.frame(cbind(1:n, log(h)))
names(df) <- c("t","lh")

m1 <- lm(formula = lh ~ t, data = df)
summary(m1)

par(opar)
plot(log(h), type = "l", main = "Regression of log stage vs. time", 
     xlab = "time index", ylab = "log stage")
abline(m1)

ld <- diff(log(h), differences = 1) 

plot(ld, type = "l", main = "differenced log stage")

# Differenced vs. detrended vs. original
par(mfrow = c(3,1))
plot(t, log(h), type = "l", main = "Original")
plot(t, resid(m1), type = "l", main = "De-trended")
plot(t[1:n-1], ld, type = "l", main = "Differenced")

# Lagged scatterplots for the monthly log data
library(astsa)
lag1.plot(log(h), max.lag = 9, corr = TRUE)

# PACF, ACF
par(mfrow = c(2,1))
acf(log(h), main = "ACF of monthly log q")
pacf(log(h), main = "PACF")

# ------------------------------------------------------------------------------------
# De-season the data and plot the acf, pacf, and spectrum

# Compute mean and standard deviation for each month
y <- as.integer(substr(as.character(t), 1, 4))
m <- as.integer(substr(as.character(t), 6, 7))
d <- as.integer(substr(as.character(t), 9, 10))
df <- as.data.frame(cbind(y, m, d, h = log(h)))

monthly.mean <- vector(length = 12)
monthly.sd <- vector(length = 12)
for (i in 1:12)
{
  monthly.mean[i] <- mean(df$h[which(df$m == i)])
  monthly.sd[i] <- sd(df$h[which(df$m == i)]) # sd follows a similar trend to the mean
}

# Subtract the mean for each month
lq.ds <- vector(length = length(df$h)) # monthly flow data, after taking log and de-seasoning
for (ind in 1:length(df$h))
{
  msubt <- monthly.mean[df$m[ind]] # choose the monthly mean to subtract
  sdsubt <- monthly.sd[df$m[ind]]
  # lq.ds[ind] <- l.a.monthly[ind] - msubt
  lq.ds[ind] <- (df$h[ind] - msubt)/sdsubt
}
df$m[ind]

# Plot monthly mean and sd (following https://www.r-bloggers.com/r-single-plot-with-two-different-y-axes/)
d = data.frame(mm =seq(1,12),
               mean = monthly.mean,
               sd = monthly.sd)
par(opar)
with(d, plot(mm, mean, type="h", col="red3", 
             ylab="Mean"),
     ylim=c(1,6), main = "Monthly Averages")
par(new = T)
with(d, plot(mm, sd, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2),
     main = "Monthly Averages")
axis(side = 4)
mtext(side = 4, line = 3, 'Standard deviation')
legend("bottomleft",
       legend=c("Mean", "SD"),
       lty=c(1,0), pch=c(NA, 16), col=c("red3", "black"))

par(mfrow = c(2,1))
plot(log(h), type = "l", main = "log stage")
plot(lq.ds, type = "l", main = "De-seasoned log stage")

mean(lq.ds) # mean is 0
sd(lq.ds) # standardized to 1

# ------------------------------------------------------------------------------------------------------------

# Compare it to a standard normal distribution
 
plot_normal <- function(x)
{
  # plots the empirical distribution and 
  # compares it to the normal distribution
  plot(density(x), type = "l", lwd = 2, ann = FALSE)
  hist(x, prob = TRUE, add = TRUE)
  t1 = seq(min(x),max(x),length=1000)
  lines(t1, dnorm(t1, mean = 0, sd = 1), type="l", col = "blue")
}

par(opar)
plot_normal(lq.ds)
legend(locator(1), c("Empirical","Theoretical N(0,1)"), lty = c(1,1), 
       cex = 0.8, col = c("black", "blue"))

plot_normal(diff(lq.ds))
plot_normal(ld)

library(MASS)
plot_gamma <- function(x)
{
  pars <- fitdistr(x, "lognormal")
  plot(density(x), type = "l", lwd = 2, ann = FALSE)
  hist(x, prob = TRUE, add = TRUE)
  t1 = seq(min(x),max(x),length=1000)
  lines(t1, dnorm(t1, mean = 0, sd = 1), type="l", col = "blue")
}

# Compare to gamma distribution
lh <- log(h)

# Try skew normal (gamma requires positive values)

rescale01 <- function(x)
{
  # rescale to 0-1
  (x-min(x))/(max(x)-min(x))
}

z <- rescale01(lh)
z[which(z==0)] <- 1e-3 # remove zero values for the fit

hist(z, breaks = "fd")
pars <- fitdistr(z, "gamma")
M <- 1000
sim.vals <- dgamma(seq(0, 10, length.out = 1000), 
                   shape = as.numeric(pars$estimate[1]), 
                   rate = as.numeric(pars$estimate[2]))
hist(sim.vals, breaks = 100, xlim = c(0,1))

# Standardize first
z <- (lh-mean(lh))/sd(lh)
hist(z,"fd")
z.rs <- rescale01(z)
hist(z.rs, "fd")

# Logistic distribution: symmetric, with fatter tails than normal.
# May be a good model for diff(lq.ds)

z <- diff(lq.ds)
pars <- fitdistr(z, "logistic")

plot(density(z), col = "red", lwd = "2", main = "Miss. River at Quincy")
hist(z, add = TRUE, breaks = 50, prob = TRUE)
distvals <- dlogis(seq(-5,5, length.out = 1e3), 
                   location = as.numeric(pars$estimate[1]), 
                   scale = as.numeric(pars$estimate[2])
                   )
lines(seq(-5,5, length.out = 1e3), distvals, col = "blue", lwd = 2)
legend("topright", 
       legend = c("Empirical distribution","Logistic distribution"), 
       fill = c("red","blue")
       )

plot(lq.ds, type = "l")
plot(z, type = "l") # one difference
plot(diff(z), type = "l") # two differences

hist(lq.ds, "fd")
hist(z, "fd")
hist(diff(z), "fd")


# -----------------------------------------------------------------------
# De-trending and de-seasonalizing

# Using detrending

y <- residuals(m1) # the residuals from the linear fit to the log stage data

# Deseasonalize it
# Compute mean and standard deviation for each month
yr <- as.integer(substr(as.character(t), 1, 4))
m <- as.integer(substr(as.character(t), 6, 7))
d <- as.integer(substr(as.character(t), 9, 10))
df <- as.data.frame(cbind(yr, m, d, h = y))

monthly.mean <- vector(length = 12)
monthly.sd <- vector(length = 12)
for (i in 1:12)
{
  monthly.mean[i] <- mean(df$h[which(df$m == i)])
  monthly.sd[i] <- sd(df$h[which(df$m == i)]) # sd follows a similar trend to the mean
}

# Subtract the mean for each month
lq.ds <- vector(length = length(df$h)) # monthly flow data, after taking log and de-seasoning
for (ind in 1:length(df$h))
{
  msubt <- monthly.mean[df$m[ind]] # choose the monthly mean to subtract
  sdsubt <- monthly.sd[df$m[ind]]
  # lq.ds[ind] <- l.a.monthly[ind] - msubt
  lq.ds[ind] <- (df$h[ind] - msubt)/sdsubt
}

# -----------------------------------------------------------------------

par(opar)
par(mfrow = c(2,1))
acf(lq.ds, lag.max = 100)
pacf(lq.ds, lag.max = 100)

acf(diff(lq.ds)) # differencing seems like it would make the residuals look like white noise.
pacf(diff(lq.ds))
# However, I did not find a significant trend in the log flow plot.

lq.dds <- diff(lq.ds) # take a difference because the log deseasoned q still doesn't look stationary

par(mfrow = c(2,1), mar = c(2,2,2,2))
plot(lq.ds, type = "l", main = "De-seasoned log q")
plot(lq.dds, type = "l", main = "De-seasoned, differenced log q")

# %%%%%%%%%%%%%%%%%%%%

# Raw periodogram
periodogram <- mvspec(lq.ds, log = "no", main = "Raw periodogram")

# Parametric periodogram fit
spaic <- spec.ar(lq.ds, log = "no")

# Raw periodogram
periodogram <- mvspec(lq.dds, log = "no", main = "Raw periodogram (differenced")
# Looks a bit like white noise to me.

# Parametric periodogram fit
spaic <- spec.ar(lq.dds, log = "no", main = "differenced")

##############################################
##############################################
# ARIMA Modeling
##############################################
##############################################

# ------------------------------------------------------
# Fit ARIMA models

# Do for each and compare results
z <- lq.ds # this distribution is much more normal than others so far
# Distributino of diff(z) is symmetric, but with a much higher peak than normal.

# Most likely models based on my assessment of ACF, PACF, spectrum
m110 <- sarima(z, 1,1,0)

m <-sarima(z, 15, 2, 0)

m1 <-sarima(z, 1, 1, 1)
m2 <-sarima(z, 1, 1, 2)
m3 <-sarima(z, 1, 1, 3)
m314 <-sarima(z, 3, 1, 4) # seems to be best

m <- arima(z, c(10,2,0))

# Compare AIC values
m1$AIC
m2$AIC
m3$AIC
m314$AIC
m110$AIC

# Compare BIC values
m1$BIC
m2$BIC
m3$BIC
m314$BIC
m110$BIC

# The ARIMA(3,1,4) model performs best

# ---------------------------------------------------------------------
# Predict stage values using the model - backtransform

set.seed(704753262)
N <- 10*365 # simulate 10 years of daily data
newdata <- arima.sim(model = mm1, n = N)
# newdata.mat <- array(data = newdata, dim = c(10,365)) # each row is one year


h.pred <- vector(length = N)
month.ind <- 1 # starts at January 1
day.ind <- 1  

for (i in 1:N)
{
  h.pred [i] <- newdata[i]*monthly.sd[month.ind]+monthly.mean[month.ind]
  
  day.ind <- day.ind + 1
  if (month.ind == 2 & day.ind > 28)
  {
    month.ind <- month.ind + 1 # ,6,9,11
  } else if ((month.ind == 4 & day.ind > 30) | 
             (month.ind == 6 & day.ind > 28) | (month.ind == 9 & day.ind > 28)
             | (month.ind == 11 & day.ind > 28))
  {
    month.ind <- month.ind + 1
  } else if (day.ind > 31)
  {
    month.ind <- month.ind + 1
  }
  if (month.ind > 12) 
  {
    month.ind <- 1
  }
}

plot(h.pred, type = "l")

# Add the trend back in
a <- coef(m1)[1]
b <- coef(m1)[2]
t2 <- (1:N) # not quite right, perhaps?
hl <- h.pred + a + b*t2
H.pred <- exp(hl)

summary(H.pred)
summary(h)

# Decent model for the first through third quantiles, but does not do
# a good job at the extreme ends of the distribution. 
# Too many low/high extreme values.

plot(h[1:1000], main = "original data", type = "l")
plot(H.pred[1:1000], main = "simulated data", type = "l")

# Does not really capture the autocorrelation present in the original data.
# Need a much higher order AR model, perhaps. AR(90) for seasons?

# Histograms
hist(h, main = "true histogram", breaks = "fd", xlim = c(138, 148))
hist(H.pred, main = "simulated histogram", breaks = "fd", xlim = c(138, 148))
# It's actually pretty good if it were truncated at min(h)
# Let's try that. But use the fractional WSE/WSE_max values, not WSE alone.

# Seasonal ARIMA model?
# High order AR model?

m.high <- arima(z, c(30,1,0))
AIC(m.high)

mm <- sarima(z,1,1,1,1,1,1,12)

acf(z)
str(pacf(z))
str(acf(diff(z, differences = 2)))
str(pacf(diff(z, differences = 2)))


mm1 <- sarima(z, 1, 1, 4)
mm2 <- sarima(z, 3, 2, 4)
mm3 <- sarima(z, 3, 1, 4)
mm1$BIC
mm2$BIC
mm3$BIC

# Try the mm1 model.

# ------------------------------------------------------------
# Long memory time series/fractional differencing

library(fracdiff)
z.fd = fracdiff(z, nar=1, nma=4, M=30)
z.fd$d # = 0.5
z.fd$stderror.dpq # = 1.79e-5 (questionable result!!)
p = rep(1,31)
for (k in 1:30){ p[k+1] = (k-z.fd$d)*p[k]/(k+1) }
plot(1:30, p[-1], ylab=expression(pi(d)), xlab="Index", type="h")
res.fd = diffseries(z, z.fd$d) # frac diff resids
res.arima = resid(arima(z, order=c(1,1,4))) # arima resids
par(mfrow=c(2,1))
acf(res.arima, 100, xlim=c(4,97), ylim=c(-.2,.2), main="")
acf(res.fd, 100, xlim=c(4,97), ylim=c(-.2,.2), main="")
# From TSA textbook, but does not give good results here.


