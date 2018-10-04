nlfit1 <- readRDS("/Users/jschap/Desktop/Cross_Sections/Outputs/pool_21_ra_TRUE_nr_3774_expo_19_spacing_5_sampling_even_mc_replicates_100/nl/nl_r_1.rds")
nlfit2 <- readRDS("/Users/jschap/Desktop/Cross_Sections/Outputs/pool_21_ra_TRUE_nr_3774_expo_19_spacing_5_sampling_even_mc_replicates_100/nl/nl_r_1_test.rds")

k <- 19
m <- 1

nlfit1[[k]][[m]] # this is much smaller
nlfit2[[k]][[m]] # this is much larger

# parameter values differ quite a bit with initial guesses

nlfit1[[k]][[m]]$call
nlfit2[[k]][[m]]$call

# Why are they different sizes?

nacount1 <- 0
nacount2 <- 0
for (k in 1:n_exp_levels)
{
  for (m in 1:M)
  {
    if(is.na(nlfit1[[k]][[m]]))
    {
      nacount1 <- nacount1 + 1
    }
    if(is.na(nlfit2[[k]][[m]]))
    {
      nacount2 <- nacount2 + 1
    }
  }
}

