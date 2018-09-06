# Monte Carlo experiments
# Takes the place of lines 167-295 in main_d4
#
# The code developed here has been implemented in main_d5 9/6/2018, though there may still be some bugs to work out.

# -----------------------------------------------------------------------------------------------------
# Do for one cross section, for now
nr <- 1
M <- 1e3 # number of replicates

expo <- seq(0.05, 0.95, length.out = n_exp_levels) # exposure levels
n_exp_levels <- length(expo)

# Initialize outputs
WSEw_obs <- vector(length = nr, "list")
lf <- vector(length = nr, "list")
for (r in 1:nr)
{
  WSEw_obs[[r]] <- vector(length = n_exp_levels, "list")
  lf[[r]] <- vector(length = n_exp_levels, "list")
}
for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    WSEw_obs[[r]][[k]] <- vector(length = M, "list")
    lf[[r]][[k]] <- vector(length = M, "list")
  }
}

begin.time <- Sys.time()
for (r in 1:nr) # loop over reaches
{
  for (k in 1:n_exp_levels) # loop over exposure levels
  {
    for (m in 1:M)
    {
      WSEw_obs[[r]][[k]][[m]] <- observe(WSEw = rWSEw[[r]], exposure = expo[k], sd_w = 0)
    }
    print(k)
  }
  if (r%%5 == 0)
  {
    current.time <- Sys.time()
    te <- current.time - begin.time
    print(paste("Processed", r, "of", nr, "cross sections")) # display progress
    print(paste("Time elapsed:", te, "units"))
  }
}

# saveRDS(WSEw_obs, file = file.path(exp_dir, "WSEw_obs.rds"))

# Fit linear model
for (r in 1:nr)
{
  for (k in 1:n_exp_levels)
  {
    for (m in 1:M)
    {
      lf[[r]][[k]][[m]] <- fit_linear(WSEw_obs[[r]][[k]][[m]])
    }
    print(k)
  }
}

# saveRDS(lf, file.path(exp_dir, "lf.rds"))

# --------------------------------------------------------------------------------------------------------
# Compute hydraulic parameters

z0.l <- array(dim = c(nr, n_exp_levels, M))
#A.l <- array(dim = c(nr, n_exp_levels, M))
#WP.l <- array(dim = c(nr, n_exp_levels, M))
A0.l <- array(dim = c(nr, n_exp_levels, M))
# Linear
for (r in 1:nr)
{
  for (k in 3:n_exp_levels)
  { 
    for (m in 1:M)
    {
      if (class(lf[[r]][[k]][[m]]) == "lm")
      {
        z0.l[r,k,m] <- predict(lf[[r]][[k]][[m]], newdata = data.frame(w = 0))
        #A.l[r,k,m] <- calc_model_A(lf[[r]][[k]][[m]], type = "linear")
        #WP.l[r,k,m] <- calc_model_WP(lf[[r]][[k]][[m]], type = "linear")
        A0.l[r,k,m] <- calc_model_A0(lf[[r]][[k]][[m]], type = "linear")
      }
    }
    print(k)
  }
  if (r %% 10 == 0)
  {
    print(paste("progress:", r, "of", nr))
  }
}
s0.l <- apply(z0.l, c(2,3), diff)

# --------------------------------------------------------------------------------------------------------
# Check how distribution of results changes with sample size

z5 <- z0.l[1,k,] # 1000

par(mfrow=c(2,2))
k <- 16
hist(z1, prob = TRUE, main = "z0, 10000 replicates, 80% exposure", xlim = c(132, 135), ylim = c(0,3), "fd")
hist(z2, prob = TRUE, main = "z0, 1000 replicates, 80% exposure", xlim = c(132, 135), ylim = c(0,3), "fd")
hist(z3, prob = TRUE, main = "z0, 500 replicates, 80% exposure", xlim = c(132, 135), ylim = c(0,3), "fd")
hist(z4, prob = TRUE, main = "z0, 100 replicates, 80% exposure", xlim = c(132, 135), ylim = c(0,3), "fd")
hist(z5, prob = TRUE, main = "z0, 50 replicates, 80% exposure", xlim = c(132, 135), ylim = c(0,3), "fd")

A5 <- A0.l[1,k,] # 1000

hist(A1, prob = TRUE, main = "A0, 10000 replicates, 80% exposure", xlim = c(100, 400), ylim = c(0,0.025), "fd")
hist(A2, prob = TRUE, main = "A0, 1000 replicates, 80% exposure", xlim = c(100, 400), ylim = c(0,0.025), "fd")
hist(A3, prob = TRUE, main = "A0, 500 replicates, 80% exposure", xlim = c(100, 400), ylim = c(0,0.025), "fd")
hist(A4, prob = TRUE, main = "A0, 100 replicates, 80% exposure", xlim = c(100, 400), ylim = c(0,0.025), "fd")
hist(A5, prob = TRUE, main = "A0, 50 replicates, 80% exposure", xlim = c(100, 400), ylim = c(0,0.025), "fd")

# --------------------------------------------------------------------------------------------------------
# Check theoretical vs. MC distribution for special case of no w error

# Theoretical distribution of z0 without measurement error
WSEw_obs1 <- observe(WSEw = rWSEw[[1]], exposure = expo[k], sd_wse = 0.1, sd_w = 0)
lf1 <- fit_linear(WSEw_obs1)
# mu <- min(rWSEw[[1]]$WSE) # unbiased assumption is incorrect because we're extrapolating...
mu <- as.numeric(lf1$coefficients[1])
sigma <- (0.1^2 + vcov(lf1)[1,1])^0.5

x <- seq(mu-6*sigma, mu+6*sigma, length=100)
hx <- dnorm(x, mean = mu, sd = sigma)
plot(x, hx, ylab="Density", col = "red", ylim = c(0,5), type = "l", main = "Hist of predicted z0 and theoretical z0 distribution")

# Histogram of predicted z0
k <- 16 # 80% exposure
hist(z0.l[1,k,], prob = TRUE, add=TRUE)
