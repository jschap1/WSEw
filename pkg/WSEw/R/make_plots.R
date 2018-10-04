# Makes figures for bathymetry project
# 
# Created 9/27/2018 JRS

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Figure 1

#' Sample Size
#' 
#' Tests appropriateness of sample size for Monte Carlo simulation
#' @details Creates histograms of z0 and A0 estimates vs. sample size for several exposure levels
#' @export
#' @param x array containing predictions of a quantity, such as z0, for a cross section
#' @examples 
#' exp_dir <- "/Volumes/HD3/Cross_Sections/pool_21_ra_10km_nr_3_spacing_5_sampling_even_MC_replicates_500"
#' n_exp_levels <- 19
#' expo <- seq(0.05, 0.95, length.out = n_exp_levels)
#' MM <- c(5,10,50,100) # number of replicates
#' r <- 1
#' k <- 4
#' z0.l <- vector(length = 4, "list")
#' par(mfrow = c(2,2))
#' for (i in 1:4)
#' {
#'   M <- MM[i]
#'   pred_lf <- pred_linear_par(r) # assumes you've already got the fitted models with the appropriate number of replicates
#'   z0.l[[i]] <- pred_lf$z0
#'   sample_size_plot(z0.l[[i]][k,], 
#'                    main = paste("z0", M, "replicates", 100*expo[k], "percent exposure"),"fd")
#' }

sample_size_plot <- function(x, ...)
{
  n <- length(x) # number of MC replicates
  hist(x, prob = TRUE, ...)
  return(0)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Figure 2

#' Monte Carlo convergence theory
#' 
#' Tests whether the Monte Carlo simulations converge to the 
#' true distribution for the linear case with no measurement error
#' @export
#' @details It is not clear whether this procedure is correct. More details should be worked out.
#' Also, this code does not work right now, but it shouldn't be too hard to get running. 
#' Just not a priority.

# Theoretical distribution of z0 without measurement error
# exp_dir <- "/Volumes/HD3/Cross_Sections/pool_21_ra_10km_nr_3_spacing_5_sampling_even_MC_replicates_500"
# load(file.path(exp_dir, "processed_xs_data_10km.rda"))
# n_exp_levels <- 19
# expo <- seq(0.05, 0.95, length.out = n_exp_levels)
# k <- 16
# r <- 1
# 
# WSEw_obs1 <- observe(WSEw = rWSEw[[1]], exposure = expo[k], sd_wse = 0.1, sd_w = 0)
# lf1 <- fit_linear(WSEw_obs1)
# 
# mu <- as.numeric(lf1$coefficients[1])
# sigma <- (0.1^2 + vcov(lf1)[1,1])^0.5
# 
# x <- seq(mu-6*sigma, mu+6*sigma, length=100)
# hx <- dnorm(x, mean = mu, sd = sigma)
# plot(x, hx, ylab="Density", col = "red", ylim = c(0,5), type = "l", main = "Hist of predicted z0 and theoretical z0 distribution")
# 
# # Histogram of predicted z0
# k <- 16 # 80% exposure
# hist(z0.l[1,k,], prob = TRUE, add=TRUE)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Figure 3

#' Fit models to fully-known raech-averaged channel geometry
#' 
#' Fits WSE-w models given fully exposed cross section geometry with no measurement error
#' @export
#' @param savename filename to save the fitted models, goodness of fit metrics, and shape parameters
#' @details
#' fits h-w models to each cross section, assuming full exposure and no measurement error
#' computes goodness of fit statistics for each model type, across all cross sections, and plots them in histograms
#' computes shape parameters for each cross section and plots them in a histogram
#' saves numerical results
#' @example fit_full_reach("reach_chars.rda")

fit_full_channel <- function(savename, plotflag = FALSE)
{
  
  nr <- length(rWSEw)
  lf <- vector(length = nr, "list")
  sb <- vector(length = nr, "list")
  sbm <- vector(length = nr, "list")
  nl <- vector(length = nr, "list")
  nlsb <- vector(length = nr, "list")
  s <- vector(length = nr) # shape parameter
  
  for (r in 1:nr)
  {
    WSEw_obs <- observe(rWSEw[[r]], exposure = 1, sd_wse = 0, sd_w = 0)
    lf[[r]] <- fit_linear(WSEw_obs)
    sb[[r]] <- fit_slopebreak(WSEw_obs, multiple_breaks = FALSE, continuity = TRUE)
    sbm[[r]] <- fit_slopebreak(WSEw_obs, multiple_breaks = TRUE, continuity = TRUE)
    nl[[r]] <- fit_nonlinear(WSEw_obs)
    nlsb[[r]] <- fit_nlsb(WSEw_obs)
    
    s[r] <- as.numeric(coef(nl[[r]])[3])
    
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
  
  if (plotflag == TRUE)
  {
    png("gof_lf.png")
    plot_gof(gof.lf)
    dev.off()
    
    png("gof_sb.png")
    plot_gof(gof.sb)
    dev.off()
    
    png("gof_sbm.png")
    plot_gof(gof.sbm)
    dev.off()
    
    png("gof_nl.png")
    plot_gof(gof.nl)
    dev.off()
    
    png("gof_nlsb.png")
    plot_gof(gof.nlsb)
    dev.off()
    
    png("shape_pars.png")
    hist(s, main = "Reach-averaged cross sections", ...)
    dev.off()
    
    print("Saved figures:")
    print("gof_lf.png")
    print("gof_sb.png")
    print("gof_sbm.png")
    print("gof_nl.png")
    print("gof_nlsb.png")
    
  }

  save(s, lf, sb, sbm, nl, nlsb, gof.lf, gof.sb, gof.sbm, gof.nl, gof.nlsb, 
       file = savename)

  print(paste("Saved fitted models, goodness of fit metrics, and shape parameters to", savename))

  gof <- list(gof.lf, gof.sb, gof.sbm, gof.nl, gof.nlsb)
  
  return(gof)
}

# ----------------------------------------------------------------------

#' Plot goodness-of-fit metrics for a particular class of model
#' 
#' Used by fit_full_channel
#' @param gof output of gof_batch
#' @param ... other arguments as taken by hist
plot_gof <- function(gof, ...)
{
  summary(gof)
  par(mfrow = c(2,3))
  hist(gof$SSE, main = "SSE", col = "darkblue", ...)
  hist(gof$AIC, main = "AIC", col = "darkblue", ...)
  hist(gof$BIC, main = "BIC", col = "darkblue", ...)
  hist(gof$r2, main = "r2", col = "darkblue", ...)
  hist(gof$MAE, main = "Mean absolute error", col = "darkblue", ...)
  hist(gof$mAE, main = "Max absolute error", col = "darkblue", ...)
}


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Figure 4

#' Characterize channel
#' 
#' Characterizes channel geometry with shape parameter and bankfull parameters
#' @export
#' @param savename filename to save the fitted models, goodness of fit metrics, and shape parameters
#' @details
#' Reports the same metrics as Grimaldi et al. (2018)
#' Drops the first d-w data point to perform fit and avoid an issue with log(0)
#' @example characterize_channel("channel_chars.rda")

characterize_channel <- function(savename, plotflag = FALSE)
{
  
  n.xs <- length(xWSEw)
  
  dist_downstream <- seq(5, n.xs*5, by = 5) # m
  wbf <- unlist(lapply(cross_sections$x, max)) # bankfull width
  dbf.max <- unlist(lapply(cross_sections$d, max)) # maximum bankfull depth
  dbf <- unlist(lapply(cross_sections$d, mean)) # average bankfull depth
  b.min <- unlist(lapply(cross_sections$b, min)) # minimum bed elevation
  
  # shape parameter (using linearized fitting method to avoid sensitivity to initial guesses)
  xdw <- vector(length = n.xs, "list")
  A <- vector(length = n.xs)
  for (r in 1:n.xs)
  {
    width <- xWSEw[[r]]$w
    depth <- xWSEw[[r]]$WSE - b.min[r]
    A[r] <- calc_A_from_WSEw(xWSEw[[r]]) # bankfull flow area
    xdw[[r]] <- data.frame(w = width/wbf[r], d = depth/dbf.max[r]) # width-depth data
  }
  
  # perform fits
  power_model <- vector(length = n.xs, "list")
  s <- vector(length = n.xs)
  for (r in 1:n.xs)
  {
    dw <- xdw[[r]][-1,] # drop the first entry to account for zero values
    power_model[[r]] <- lm(log(d) ~ log(w) + 0, data = dw)
    s[r] <- as.numeric(coef(power_model[[r]]))
  }
  
  if (plotflag == TRUE)
  {
    
    png("cross_section_pars_1.png")
    par(mfrow = c(4,1))
    plot(dist_downstream/1000, dbf, type = "l", 
         main = "maximum bankfull depth (m)", xlab = "distance downstream (km)",
         ylab = "dbf_max")
    plot(dist_downstream/1000, dbf.max, type = "l", 
         main = "maximum bankfull depth (m)", xlab = "distance downstream (km)",
         ylab = "dbf_max")
    plot(dist_downstream/1000, wbf, type = "l", 
         main = "bankfull width (m)", xlab = "distance downstream (km)",
         ylab = "wbf")
    plot(dist_downstream/1000, b.min, type = "l", 
         main = "minimum bed elevation (m)", xlab = "distance downstream (km)",
         ylab = "b_min")
    dev.off()
    
    png("cross_section_pars_2.png")
    par(mfrow = c(2,1))
    hist(s, main = "shape parameters")
    plot(dist_downstream/1000, A, type = "l", 
         main = "bankfull flow area (sq. m)", 
         xlab = "distance downstream (km)",
         ylab = "b_min")
    dev.off()
    
    print("Saved plots as cross_section_pars_1 and cross_section_pars_2.")
  }
  
  save(A, s, wbf, dbf, dbf.max, dist_downstream, 
       b.min, xdw, power_model, 
       file = savename)
  
  print("Saved cross section parameters as savename.")
  
}

# ------------------------------------------------------------------------------
# Figure 5

#' Hydraulic parameter predictions boxplots
#' 
#' Makes box and whisker plots showing predicted hydraulic parameter error
#' @param z list of arrays containing parameter predictions from each class of model
#' @param z.true true value of the parameter
#' @param k exposure level
#' @param lumped flag, if true, lump the cross sections together and just compare methods
#' @param legend flag for legend
#' @param absolute flag for plotting absolute error
#' @param ... other parameters as taken by boxplot()
#' @export
#' @details If lumped = FALSE, plots the results over each cross section separately.
#' @examples 
#' expo <- seq(0.05, 0.95, length.out = 19)
#' k <- 10
#' parameter_predictions_boxplots(z0.l, z0.sb, z0.nl, z0.nlsb, z0.true.ra, k = k, lumped = TRUE, 
#'                                main = paste("z0 error for each method at", 100*expo[k], "percent channel exposure"), 
#'                                ylab = "z0 error (m)")
#' parameter_predictions_boxplots(z0.l, z0.sb, z0.nl, z0.nlsb, z0.true.ra, k = 12, absolute = FALSE, lumped = TRUE)
#' par(mfrow = c(2,2))
#' k <- 16
#' parameter_predictions_boxplots(z0.l, z0.sb, z0.nl, z0.nlsb, z0.true.ra, k = k, lumped = FALSE, absolute = FALSE, 
#'                                main = paste("z0 error at", 100*expo[k], "percent channel exposure"), 
#'                                ylab = "z0 error (m)", legend = FALSE)  

bplab <- c("1","1","1","1","2","2","2","2","3","3","3","3")
par(mfrow = c(2,2))
k <- 16
parameter_predictions_boxplots(A0.l, A0.sb, A0.nl, A0.nlsb, A0.true.ra, k = k, lumped = FALSE, absolute = FALSE,
                               main = paste("A0 error at", 100*expo[k], "percent channel exposure"),
                               ylab = "A0 error (sq. m)", legend = FALSE, A0 = TRUE, notch = TRUE, names = bplab)

# parameter_predictions_boxplots(z0.l, z0.sb, z0.nl, z0.nlsb, z0.true.ra, k = k, lumped = FALSE, absolute = FALSE,
#                                main = paste("z0 error at", 100*expo[k], "percent channel exposure"),
#                                ylab = "z0 error (m)", legend = FALSE, A0 = FALSE, notch = TRUE, names = bplab)

parameter_predictions_boxplots <- function(z.l, z.sb, z.nl, z.nlsb, z.true, k,
                                           lumped = FALSE, absolute = FALSE,
                                           legend = TRUE, A0 = FALSE, ...)
{
 
  z <- c(as.vector(t(z.l[,k,])),
         as.vector(t(z.sb[,k,])),
         as.vector(t(z.nl[,k,])),
         as.vector(t(z.nlsb[,k,])))
  
  if (A0)
  {
    z.t <- rep(rep(z.true[,k], each = 500), 4)
  } else
  {
    z.t <- rep(rep(z.true, each = 500), 4)
  }
  
  
  z.err <- z - z.t
  
  if (absolute)
  {
    z.err <- abs(z.err)
  }
  
  xs <- rep(rep(1:3, each = 500), 4)
  
  type <- rep(c(rep("L", 500),
                rep("SB", 500),
                rep("NL", 500),
                rep("NLSB", 500)), 3)
  
  df <- data.frame(z.err = z.err, xs = xs, type = type)
  df$type=factor(df$type, levels=levels(df$type)[c(1,4,2,3)]) # reorder
  
  cols <- c("red","orange","green","blue")
  
  if (lumped)
  {
    boxplot(z.err ~ type, df, col = cols, ...)
    abline(0,0, lty = 2, lwd = 1.5)

  } else
  {
    boxplot(z.err ~ xs + type, df, col = cols, ...)
    abline(0,0, lty = 2, lwd = 1.5)
    if (legend)
    {
      legend("topleft", horiz = TRUE, legend = c("L","SB","NL","NLSB"), fill = cols)
    }
  }

  return(0)
}

# ------------------------------------------------------------------------------
# Figure 6

#' Parameter predictions scatterplots (without error)
#' 
#' Plot predicted z0 and A0 from each model for each reach-average cross section
#' Assumes there is M = 1 replicate. This should be the no error case.
#' This is an earlier iteration of Figure 5
#' @param z0.l
#' @param z0.sb
#' @param z0.nl
#' @param z0.nlsb
#' @export

parameter_predictions_scatterplots <- function(z0.l, z0.sb, z0.nl, z0.nlsb)
{
  par(mfrow = c(2,2))
  
  for (k in seq(4,16,by=4))
  {
    plot(xlim = c(0,4), 
         ylim = c(100, 144), xlab = "cross section", ylab = "z0 (m)",
         main = paste("z0 predictions,", 100*expo[k], "percent exposure"), 
         "n")
    if (k==16)
    {
      legend("bottomright", legend = c("Truth", "Linear","SB","NL","NLSB"),
             col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)
    }
    for (r in 1:nr)
    {
      rx <- seq(r-0.25, r+0.25, by = 0.01)
      ry <- rep(z0.true.ra[r], length(rx))
      lines(rx, ry, col = "black", pch = 7)
      points(r, z0.l[r,k,1], col = "red", pch = 1, cex = 1.5)
      points(r, z0.sb[r,k,1], col = "orange", pch = 1, cex = 1.5)
      points(r, z0.nl[r,k,1], col = "green", pch = 1, cex = 1.5)
      points(r, z0.nlsb[r,k,1], col = "blue", pch = 1, cex = 1.5)
    }
  }
  
  par(mfrow = c(2,2))
  for (k in seq(4,16,by=4))
  {
    plot(xlim = c(0,4), 
         ylim = c(0, 4500), xlab = "cross section", ylab = "A0 (sq. m)",
         main = paste("A0 predictions,", 100*expo[k], "percent exposure"), 
         "n")
    if (k==16)
    {
      legend("topright", legend = c("Truth", "Linear","SB","NL","NLSB"),
             col = c("black", "red","orange","green","blue"), lwd = c(1,1,1,1,1), ncol = 3)
    }
    for (r in 1:nr)
    {
      rx <- seq(r-0.25, r+0.25, by = 0.01)
      ry <- rep(A0.true.ra[r,k], length(rx))
      lines(rx, ry, col = "black", pch = 7)
      points(r, A0.l[r,k,1], col = "red", pch = 1, cex = 1.5)
      points(r, A0.sb[r,k,1], col = "orange", pch = 1, cex = 1.5)
      points(r, A0.nl[r,k,1], col = "green", pch = 1, cex = 1.5)
      points(r, A0.nlsb[r,k,1], col = "blue", pch = 1, cex = 1.5)
    }
  }
  return(0)
}



# ------------------------------------------------------------------------------
# Figure 7



# ------------------------------------------------------------------------------
# Figure 7



# ------------------------------------------------------------------------------
# Figure 7



# ------------------------------------------------------------------------------
# Figure 7






# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SCRAP

# length(rep(rep(1:3, each = 500), 4))

# z1 <- t(z0.l[,k,])
# summary(z1[,2])
# z2 <- as.vector(z1)
# summary(z2[501:1000]) # this demonstrates how as.vector works

# z0 <- c(as.vector(t(z0.l[,k,])),
#   as.vector(t(z0.sb[,k,])),
#   as.vector(t(z0.nl[,k,])),
#   as.vector(t(z0.nlsb[,k,])))

# cols <- rainbow(4, s = 0.5)

# by reach

