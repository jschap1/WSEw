# Executive file for WSEw project
# 
# Created 8/2/2018 JRS
# Revised 8/10/2018 JRS
# Revised 9/5/2018 JRS 
#   Reorganized for ease of use
#   Removed scrap code from end of document
#   Changed looping structure to avoid having to load all models into memory at once
# Revised 9/6/2018 JRS
#   Added Monte Carlo simulation capability. Can set M = 1 replicate for deterministic case.
#   Changed the way models are stored to avoid having to store too much data in RAM.
# Revised 9/13/2018
#   Implemented multicore processing with foreach package

# Parallel processing -------------------------------------------------------------------------------------------
library(foreach)
library(doMC)
library(doRNG)

registerDoMC(cores = 2)
seed <- doRNGseed() # this function cannot be found for some reason...

# Example for reproducibility
a <- foreach(i=1:2, .combine = cbind) %dorng% {rnorm(5)}
doRNGseed(seed)
b <- foreach(i=1:2, .combine = cbind) %dorng% {rnorm(5)}
identical(a,b)

# a <- foreach(i=1:2, .combine = cbind) %dopar% {rnorm(5)}

# -----------------------------------------------------------------------------------------------------------------------------
# Make observations

for (r in 1:nr) # loop over reaches
{
}

observe_batch <- function()
{
  WSEw_obs <- array()
  for (k in 1:n_exp_levels) # loop over exposure levels
  {
    for (m in 1:M)
    {
      WSEw_obs[[r]][[k]][[m]] <- observe(WSEw = rWSEw[[r]], exposure = expo[k])
    }
  }
  return(WSEw_obs)
}
# Save
saveRDS(WSEw_obs, file = file.path(exp_dir, "WSEw_obs.rds"))

# Leave observe() as serial. It will be easier to parallelize the function fitting because the format [[k]][[m]] is good for it.

# -----------------------------------------------------------------------------------------------------------------------------
# Fit linear model

fit_linear_par <- function(r)
{
  lf <- vector(length = n_exp_levels, "list")
  for (k in 1:n_exp_levels)
  {
    lf[[k]] <- vector(length = M, "list")
    for (m in 1:M)
    {
      model1 <- fit_linear(WSEw_obs[[r]][[k]][[m]])
      if (!is.null(model1))
      {
        lf[[k]][[m]] <- model1
      } else
      {
        lf[[k]][[m]] <- model1 <- NA
      }
    }
  }
  lf_name <- paste0("lf/lf_", "r_", r, ".rds")
  saveRDS(lf, file = file.path(exp_dir, lf_name))
  return(lf)
}

# inputs: n_exp_levels, M, r, exp_dir
nr <- 50
registerDoMC(cores = 7)
begin.time <- Sys.time()
lf1 <- foreach(r = 1:nr) %dopar% fit_linear_par(r)
print(Sys.time() - begin.time)
# There is a major benefit in terms of processing time. From 3 minutes to 40 seconds for nr=50, M = 100, fit_linear.


# Also should work out how to do this for observe() and for prediction, but the effort may not be worth it; it will only save about 2 days of computation time.



  
  






