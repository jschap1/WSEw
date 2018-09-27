# ------------------------------------------------------------------------------------------------
# Sample heights according to gauge data

# load gauge data
load(file.path("/Users/jschap/Desktop/Cross_Sections/Data", "quincy_gauge_data.rda"))
# gauge$Height # about 70 years of data, choose high, medium, and low periods
hma <- ra(gauge$Height, 3*365)

ncores <- detectCores()
registerDoMC(cores = ncores - 1)
brkpts <- breakpoints(Height~1, data = gauge, hpc = "foreach", h = 5000) # takes a long time because DP
# see http://r.789695.n4.nabble.com/Strucchange-Breakpoint-slow-td4634602.html for explanation

bp <- 12844
nn <- length(gauge$Height)
h1 <- gauge$Height[1:bp]
h2 <- gauge$Height[(bp+1):nn]
m1 <- lm(Height~1, gauge[1:bp,])
m2 <- lm(Height~1, gauge[(bp+1):nn,])

par(opar)
plot(timevector, gauge$Height, type="l")
lines(timevector[1:bp], predict(m1, newdata = data.frame(1:bp)), col = "blue")
lines(x = timevector[(bp+1):nn], y = predict(m2, newdata = data.frame((bp+1):nn)), col = "blue")

mean(h1, na.rm = TRUE) # the horizontal fits are just the mean of the data
mean(h2, na.rm = TRUE)

# Now use this approach to find some representive 3-year periods for SWOT experiment
# But breakpoints is too computationally intensive. Try a faster algorithm.
# One option is to use a coarser time grid
# PELT in the changepoints package might work. 
# There are efficient implementations of the DP algorithm for the simple location model y~1

ht <- ts(gauge$Height) # try smoothing it to get a better idea of which years are high and which are low.

# aggregate to monthly data
library(lubridate)
library(plyr) 
df <- data.frame(date = timevector, x = gauge$Height)
df$my <- floor_date(df$date, "month")
df.mon <- ddply(df, "my", summarise, x = mean(x))
plot(x~my, df.mon, type = "l")

brkpts <- breakpoints(x~1, data = df.mon, hpc = "foreach", h = 12*3) # takes about 1 minute
# One breakpoint was found at observation number 293 (months).
# The first 24 years or so are relatively low flow, and the rest are relatively high flow

# aggregate to yearly data
df <- data.frame(date = timevector, x = gauge$Height)
df$my <- floor_date(df$date, "year")
df.yr <- ddply(df, "my", summarize, x = mean(x, na.rm = TRUE))

par(mfrow = c(2,1))
plot(x~my, df.yr, type = "l")
# plot(ra(df.yr$x, 3), type=  "l")

# Based on looking at the annual average plot:
# Low 1954-1956
# Medium 1971-1973
# High 2009-2011

# plot(x~my,as.ts(df.yr))

# ------------------------------------------------------------------------------------------------
# Sample heights using the gauge data: 
# Come in with water heights and extract the widths. Do for each set of h-w data from bathymetry

startyr <- 1954
t1 <- timevector[year(timevector) >= startyr & year(timevector) <= startyr+2]
h1 <- gauge$Height[year(timevector) >= startyr & year(timevector) <= startyr+2]
plot(t1, h1, type = "l")

# fraction of bankfull height/depth at the gauge location
h.bf <- 458.59+17 # based on (incomplete) info from USACE
# throw out any out-of-bank observations for now
summary(gauge$Height/h.bf)
summary(h1/h.bf)



gauge_data <- data.frame(t = t1, h = h1)
xWSEw <- generate_hw(cross_sections, gauge_data, dx = 1)

# Method 1: find corresponding flow width for WSE ranging from empty to bankfull conditions
xWSEw <- calc_WSEw(cross_sections, interval = 0.05, dx = 1) # number of data points depends on discretization

rWSEw <- reach_avg(xWSEw)
rWSEw_burr <- reach_avg(xWSEw)

save(cross_sections, xWSEw, rWSEw, file = file.path(saveloc, "/Data/Processed_Data/processed_xs_data.rda"))

# Plot the observations
plot(WSE~w, rWSEw[[1]], main = "WSE-w sampling, three years")
points(WSE~w, rWSEw_burr[[1]], col="red", pch=19)
legend("topleft", legend = c("Even sampling","Burr sampling"), fill = c("black", "red"))