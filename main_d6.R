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
seed <- doRNGseed()

# Example for reproducibility
a <- foreach(i=1:2, .combine = cbind) %dorng% {rnorm(5)}
doRNGseed(seed)
b <- foreach(i=1:2, .combine = cbind) %dorng% {rnorm(5)}
identical(a,b)

# a <- foreach(i=1:2, .combine = cbind) %dopar% {rnorm(5)}


