# Time series modeling of river stage
# Aug. 29, 2018

setwd("C:/Users/Jacob/Box/Margulis_Research_Group/Jacob/UMBB/Data/")
load("stage_height.rda")

h <- h/max(h)

# ----------------------------------------------------------------------

# De-trend
n <- length(t)
df <- as.data.frame(cbind(1:n, log(h)))
names(df) <- c("t","lh")
m1 <- lm(formula = lh ~ t, data = df)
y <- residuals(m1) # the residuals from the linear fit to the log stage data

# De-season

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
  lq.ds[ind] <- (df$h[ind] - msubt)/sdsubt
}

# Fit ARIMA model
z <- lq.ds
m114 <- sarima(z, 1, 1, 4)

# Make predictions
set.seed(704753262)
N <- 1
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

# Add the trend back in
a <- coef(m1)[1]
b <- coef(m1)[2]
t2 <- (1:N) # not quite right, perhaps?
hl <- h.pred + a + b*t2
H.pred <- exp(hl)

plot(h[1:1000], main = "original data", type = "l")
plot(H.pred[1:1000], main = "simulated data", type = "l")

# Histograms
hist(h, main = "true histogram", breaks = "fd", xlim = c(0.9, 1))
hist(H.pred, main = "simulated histogram", breaks = "fd", xlim = c(0.9, 1))

# ----------------------------------------------------------------
# Do sampling in a function

sample.arima <- function(model = mm114, N = 1)
{
  
  newdata <- arima.sim(model, n = N)
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
  
  # Add the trend back in
  a <- coef(m1)[1]
  b <- coef(m1)[2]
  t2 <- (1:N) # not quite right, perhaps?
  hl <- h.pred + a + b*t2
  H.pred <- as.numeric(exp(hl))
  
  return(H.pred)
}

# ----------------------------------------------------------------
# Sample from the truncated distribution

M <- 1000 # number of samples to draw
pred.vals <- vector(length = M)
k <- 0
while (k<M)
{
  H.pred <- sample.arima()
  if(H.pred >= min(h))
  {
    k <- k + 1
    pred.vals[k] <- H.pred
  }
}

hist(h, main = "true histogram", 
     breaks = "fd", 
     xlim = c(0.9, 1)
     )
hist(pred.vals, main = "simulated histogram", 
     breaks = "fd", 
     xlim = c(0.9, 1)
     )

summary(h)
summary(pred.vals)

# This does a better job simulating the lower end of the distribution.
# The higher end is not so well-simulated.
# Try applying this to sample WSE from the Pool 21 cross sections.