# Builds test cases to check shape fitting algorithm

WSEref <- 25

# Parabolic cross section
x <- seq(0,100,length.out = 30)
y <- 0.01*(x-50)^2+0
plot(x,y)
parabola <- data.frame(x=x, b=y, d=WSEref-y)

# Triangular cross section
x <- seq(0,100,length.out = 30)
y1 <- 25 - 0.5*x[1:15]
y2 <- 0.5*(x[16:30]-50)
y <- c(y1, y2)
plot(x,y)
triangle <- data.frame(x=x, b=y, d=WSEref-y)

# Rectangular cross section
x <- seq(0,100,length.out = 30)
y <- rep(0,30)
y[1] <- 25
y[30] <- 25
plot(x,y)
rectangle <- data.frame(x=x, b=y, d=WSEref-y)

cross_sections <- list()
cross_sections$x <- list(rectangle$x, triangle$x, parabola$x)
cross_sections$b <- list(rectangle$b, triangle$b, parabola$b)
cross_sections$d <- list(WSEref - rectangle$b, WSEref - triangle$b, WSEref - parabola$b)

WSEw <- calc_WSEw(cross_sections, interval = 0.05, dx = 1)
save(cross_sections, WSEw, file = "test_cross_sections.rda")
