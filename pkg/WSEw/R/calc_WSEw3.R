#' Calculate WSE-w
#' 
#' Calculates WSE-w pairs for each cross section.
#' @export
#' @param cross_sections list of x coordinates, bed elevations, and depths for each transect
#' @param dist type of distribution for SWOT-sampled stage values.
#' @param pars vector of distribution parameters
#' @param n.obs number of SWOT observations to make. The UMRB is observed 1-2 times per 21-day cycle
#' @param dx used for computing width with trapezoidal rule
#' @details Makes a vector of WSE values from the minimum bed elevation to the bankfull WSE
#' Computes flow width for each WSE-w value.
#' Puts the estimated flow widths into the WSEw data frame
#' Update from old version: Now takes the SWOT sampling into account, instead of using evenly-spaced measurements
#' Exponential distribution was found to be the best fit to in-situ stage data for pool 21 of the Mississippi River
#' @return WSEw data frame containing WSE and flow width values
#' @examples 
#' xWSEw <- calc_WSEw(cross_sections, dist = "exponential", pars = c(0.9428, 74.06), n.obs = floor(1*365/10))
#' xWSEw <- calc_WSEw(cross_sections, dist = "gamma", pars = c(9022.2817, 9434.2034), n.obs = floor(1*365/10))

# OUTPUTS
# WSE-w pairs for each cross section: xWSEw, cross section geometry

calc_WSEw3 <- function(cross_sections, dist = "exponential", pars = c(0.9428, 74.06), dx = 1, n.obs)
{
  
  x <- cross_sections$x
  b <- cross_sections$b
  d <- cross_sections$d
  
  if (dist == "exponential")
  {
    shift <- pars[1] # Note that this may be valid only for a particular region (here, Pool 21)
    rate <- pars[2]
  } else if (dist == "gamma")
  {
    shape = pars[1]
    rate <- pars[2]
  }
  
  if (dist == "exponential")
  {
    x1 <- shift + rexp(n = n.obs, rate = rate) # exponential random variables telling the fraction of WSE_bf to use
  } else if (dist == "gamma")
  {
    x1 <- rgamma(n = n.obs, shape = shape, rate = rate)
  } else if(dist == "burr")
  {
    x1 <- sim_z(n = n.obs)
  } 
  x1[x1>1]<- 1 # overbank flows are not used
  
  # Pre-processing
  nseg <- length(b)
  dbf <- unlist(lapply(d, max, na.rm=TRUE)) # bankfull depth (m)
  b.min <- unlist(lapply(b, min)) # bottom depth
  
  xWSEw <- vector(length = nseg, "list")
  for (seg in 1:nseg)
  {
    
    # Make a vector of WSE(t) values to enter in with
    if (dbf[seg] < 0){next} # error check
    
    WSE_bf <- b.min[seg] + dbf[seg] # bankfull WSE
    # print(WSE_bf)
    
    WSE <- x1*WSE_bf
    #print(WSE)
    
    # Estimate the flow width for each WSE value
    #print(x[[seg]])
    #print(b[[seg]])
    w <- get_width(WSE, x[[seg]], b[[seg]], delx = dx)
    
    # Put WSE-w info into a data frame
    WSEw <- data.frame(WSE = WSE, w = w)
    
    xWSEw[[seg]] <- WSEw
  }
  
  # Remove null cross sections, those with no measurements (optional)
  rm.ind <- which(unlist(lapply(xWSEw, is.null)))
  if (any(rm.ind))
  {
    xWSEw <- xWSEw[-rm.ind] # this doesn't always seem to work...
  }
  
  return(xWSEw)
}

# ------------------------------------------------------------------------------------------------------------------------

#' Simulate SWOT WSE
#'
#' Simulate z values for SWOT sampling using a fitted Burr distribution
#' @export
#' @param scale
#' @param shape1
#' @param shape2
#' @param m linear model used to detrend the data
#' @param n number of replicates
#' @example z.sim <- sim_z(scale, shape1, shape2, m1, min(z), max(z), n=100)
#' @importFrom actuar burr

sim_z <- function(scale = 0.3187264, shape1 = 0.3211963, shape2 = 25.1872518, 
                  zmin = 0.9375944, zmax = 1,  n)
{
  z_rs <- rburr(n, shape1 = shape1, shape2 = shape2, scale = scale)
  z <- zmin + z_rs*(zmax - zmin)
}

