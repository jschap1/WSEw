# Check how well each type of model fits the true channel using true 
# measurements and a variety of model selection metrics

opar <- par()

library(WSEw)

setwd("/Users/jschap/Desktop/Cross_Sections") # best coding practice is to not use setwd...

saveloc <- "/Users/jschap/Desktop/Cross_Sections/Data/Processed_Data"
load(file.path(saveloc, "processed_data_p21_sl_5m_hires.rda"))

# If have already run, just load the results
load(file.path(saveloc, "model_selection_results.rda"))

# --------------------------------------------------------------------------------------------------
# Compute standard goodness-of-fit measures for each class of model:
# Linear model
# Slope break model (one slope break, continuous)
# Slope break model (multiple slope breaks, continuous)
# Nonlinear model
# NLSB (one "slope" break, continuous)

nr <- length(rWSEw)
lf <- vector(length = nr, "list")
sb <- vector(length = nr, "list")
sbm <- vector(length = nr, "list")
nl <- vector(length = nr, "list")
nlsb <- vector(length = nr, "list")

# linear model fits
for (r in 1:nr)
{
  WSEw_obs <- observe(rWSEw[[r]], exposure = 1, sd_wse = 0, sd_w = 0)
  lf[[r]] <- fit_linear(WSEw_obs)
  sb[[r]] <- fit_slopebreak(WSEw_obs, multiple_breaks = FALSE, continuity = TRUE)
  sbm[[r]] <- fit_slopebreak(WSEw_obs, multiple_breaks = TRUE, continuity = TRUE)
  nl[[r]] <- fit_nonlinear(WSEw_obs)
  nlsb[[r]] <- fit_nlsb(WSEw_obs)
  if (r%%5 == 0)
  {
    print(paste("Processed", r, "of", nr, "cross sections")) # display progress
  }
}

gof.lf <- gof_batch(lf, type = "l")
gof.sb <- gof_batch(sb, type = "sb")
gof.sbm <- gof_batch(sbm, type = "sbm")
gof.nl <- gof_batch(nl, type = "nl")
gof.nlsb <- gof_batch(nlsb, type = "nlsb")

save(lf, sb, sbm, nl, nlsb, gof.lf, gof.sb, gof.sbm, gof.nl, gof.nlsb, 
     file = file.path(saveloc, "model_selection_results.rda"))

# --------------------------------------------------------------------------------------------------
# Figure for paper
# Shape parameter for the reach averaged cross sections

opar <- par()
par(mfrow = c(2,1), mar = c(5,5,2,5))

nl.xs <- vector(length = n.xs, "list")
for (k in 1:n.xs)
{
  WSEw_obs <- observe(xWSEw[[k]], exposure = 1, sd_wse = 0, sd_w = 0)
  nl.xs[[k]] <- fit_nonlinear(WSEw_obs)
  if (k %% 10 == 0){print(k)}
}

s.xs <- vector(length = n.xs)
for (k in 1:n.xs)
{
  s.xs[k] <- as.numeric(coef(nl.xs[[k]])[3])  
}
summary(s.xs)

# Shape parameter for the cross sections

s <- vector(length = nr)
for (k in 1:nr)
{
  s[k] <- as.numeric(coef(nl[[k]])[3])  
}
summary(s)

hist(s.xs, main = "Cross sections", col = "blue", breaks = "stu")
hist(s, main = "Reach-averaged cross sections", col = "blue", breaks = "stu", xlim = c(0, 3.5))

# --------------------------------------------------------------------------------------------------
# Make plots

gof <- gof.sb
summary(gof)
par(mfrow = c(2,3))
hist(gof$SSE, main = "SSE", col = "darkblue")
hist(gof$AIC, main = "AIC", col = "darkblue")
hist(gof$BIC, main = "BIC", col = "darkblue")
hist(gof$r2, main = "r2", col = "darkblue")
hist(gof$MAE, main = "Mean absolute error", col = "darkblue")
hist(gof$mAE, main = "Max absolute error", col = "darkblue")




