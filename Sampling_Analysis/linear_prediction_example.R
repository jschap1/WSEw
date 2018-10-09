# Linear regression example

x <- rnorm(50)
y <- 2*x + 7 + rnorm(50, sd = 1)
plot(x,y, main = "Linear regression")
df <- data.frame(x=x, y=y)

lf <- lm(y~x, df)
summary(lf)
lines(x, predict(lf))

# Standard error of estimate

# Prediction (manual)
range(x)
x0 <- 3
a <- as.numeric(coef(lf)[1])
b <- as.numeric(coef(lf)[2])
a.var <- vcov(lf)[[1]]
b.var <- vcov(lf)[[2]]
pred.var <- (1 + a.var + x0^2*b.var)
upr <- a+b*x0 + 1.96*pred.var^0.5
lwr <- a+b*x0 - 1.96*pred.var^0.5

predict(lf, newdata = data.frame(x=x0), interval = "prediction") # (with R)

# The upper and lower bounds should match, but they don't quite.
# Investigate further.

