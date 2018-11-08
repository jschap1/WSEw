#' Characterize channel
#' 
#' Characterizes channel geometry with shape parameter and bankfull parameters
#' @export
#' @param cross_sections cross section geometry
#' @param xWSEw WSE-w data
#' @param savename filename to save the fitted models, goodness of fit metrics, and shape parameters
#' @details
#' Reports the same metrics as Grimaldi et al. (2018)
#' Drops the first d-w data point to perform fit and avoid an issue with log(0)
#' @example characterize_channel(cross_sections_avg, "channel_chars.rda")
#' @return list with the following components:
#' A, cross sectional area
#' s, shape parameter
#' wbf, bankfull width
#' dbf, mean bankfull depth
#' dbf.max, maximum bankfull depth
#' dist_downstream, distance downstream
#' b.min, minimum bed elevation
#' xdw, data frame containing width-depth data
#' power model, models fit to linearized cross section model d = aw^s

characterize_channel <- function(cross_sections, xWSEw, savename = NULL, plotflag = FALSE, section_length = 5, na.rm = FALSE)
{
  
  n.xs <- length(xWSEw)
  dist_downstream <- 1:n.xs*section_length
  na.ind <- get_na_ind(cross_sections)
  
  # # Handle NA values
  # if (na.rm)
  # {
  #   na.ind <- get_na_ind(cross_sections)
  #   cross_sections_sub <- subset_cross_sections(cross_sections, which(!na.ind))
  #   cross_sections <- cross_sections[-na.ind]
  #   xWSEw_sub <- subset_rWSEw(xWSEw, which(!na.ind))
  #   dist_downstream <- dist_downstream[-na.ind]
  # }
  
  wbf <- unlist(lapply(cross_sections$x, max, na.rm = na.rm)) # bankfull width
  dbf.max <- unlist(lapply(cross_sections$d, max, na.rm = na.rm)) # maximum bankfull depth
  dbf <- unlist(lapply(cross_sections$d, mean, na.rm = na.rm)) # average bankfull depth
  b.min <- unlist(lapply(cross_sections$b, min, na.rm = na.rm)) # minimum bed elevation
  
  wbf[na.ind] <- NA
  dbf.max[na.ind] <- NA
  dbf[na.ind] <- NA
  b.min[na.ind] <- NA
  
  # shape parameter (using linearized fitting method to avoid sensitivity to initial guesses)
  xdw <- vector(length = n.xs, "list")
  A <- vector(length = n.xs)
  power_model <- vector(length = n.xs, "list")
  s <- vector(length = n.xs)
  
  for (r in 1:n.xs)
  {
    
    if (na.ind[r])
    {
      next
    }
    
    if (dbf.max[r] == 0) # there is a division error when dbf.max is zero
    {
      next
    }
    
    width <- xWSEw[[r]]$w
    depth <- xWSEw[[r]]$WSE - b.min[r]
    A[r] <- calc_A_from_WSEw(xWSEw[[r]]) # bankfull flow area
    xdw[[r]] <- data.frame(w = width/wbf[r], d = depth/dbf.max[r]) # width-depth data; 
    
    # perform fits
    z.ind <- which(xdw[[r]]$w == 0)
    dw <- xdw[[r]][-z.ind,] # drop the first entry to account for zero values
    power_model[[r]] <- lm(log(d) ~ log(w) + 0, data = dw)
    s[r] <- as.numeric(coef(power_model[[r]]))
    
  }
  
  s[na.ind] <- NA
  A[na.ind] <- NA

  # dw <- xdw[[r]][-1,] # drop the first entry to account for zero values
  # if (r==176){dw <- dw[[r]][-1,]} # sometimes, there is a cross section with two zeros
  
  if (plotflag == TRUE)
  {
    
    png("bankfull_parameters.png")
    par(mfrow = c(2,3))
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
    plot(dist_downstream/1000, A, type = "l", 
         main = "bankfull flow area (sq. m)", 
         xlab = "distance downstream (km)",
         ylab = "A")
    dev.off()
    
    png("shape_parameters.png")
    par(mfrow = c(2,1))
    plot(dist_downstream/1000, s, type = "l", 
         main = "shape parameter", 
         xlab = "distance downstream (km)",
         ylab = "s")
    abline(1,0, col = "gray")
    abline(2,0, col = "gray")
    hist(s, main = "shape parameters", "fd", col = "lightgray")
    dev.off()
    
    print("Saved plots as bankfull_parameters.png and shape_parameters.png")
  }
  
  if (!is.null(savename))
  {
    save(A, s, wbf, dbf, dbf.max, dist_downstream, 
         b.min, xdw, power_model, 
         file = savename)
    print(paste("Saved cross section parameters as", savename)) 
  }

  result <- list(A = A, s = s, wbf = wbf, dbf = dbf, dbf.max = dbf.max, 
                 dist_downstream = dist_downstream, b.min = b.min, 
                 xdw = xdw, power_model = power_model)
  return(result)
  
}
