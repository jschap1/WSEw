# Hydraulic parameters comparison (true vs. fitted)
#
# For now, only comparing "linear"

rm(list=ls())
setwd("/Users/jschap/Desktop/Cross_Sections")

# Load cross section geometry and WSE-w relationships
processed_data <- "/Users/jschap/Desktop/Cross_Sections/Data/Processed_Data/processed_data_p21_500m.rda"
load(processed_data)

# Save names
saveloc <- "/Users/jschap/Desktop/Cross_Sections/Data/Fitting_Results/"
sbsavename <- "r_nr76_expo20_500m_0.05_sb.rda" # slope break
lsavename <- "r_nr76_expo20_500m_0.05_l.rda" # linear
nlsavename <- "r_nr76_expo20_500m_0.05_nl.rda" # nonlinear
nlsbsavename <- "xs_nr76_expo10_500m_0.05_nlsb.rda" # shape break

load(file.path(saveloc, lsavename))

# For a single reach-average cross section
r <- 1
WSEw <- rWSEw[[r]]

# --------------------------------------------------------------------------------------------------
# Geometry functions

flow_area_triang <- function(w, WSE, z0)
{
  A <- 0.5*w*(WSE - z0)
  return(A)
}

# --------------------------------------------------------------------------------------------------
# Cross sectional area

# Assume z0 is known
z0.true[r]
wse.bf <- max(WSEw$WSE)
w.bf <- max(WSEw$w)

# At bankfull

# True cross sectional area
A.bf.true <- flow_area_triang(w.bf, wse.bf, z0.true[r])

# Estimated cross sectional area
z0.est <- z0.bar[r,19] # takes the average z0 estimated for 100% channel exposure
A.bf.est <- flow_area_triang(w.bf, wse.bf, z0.est)

# The linear method underestimates z0, so it overestimates flow area of the triangular cross section
# But how does the triangular cross section compare to the true cross section, even when z0 is known?

# Cross sectional area of the full, irregular cross section

# 1. Calculate cross sectional area for each cross section in the reach. Use numerical method.
# 2. Average these cross sectional areas to get the reach-average cross section. See reach-averaging code.

cross.sections[[1]]
plot(cross.sections[[r]]$x, cross.sections[[r]]$b)

# Riemann sum (trapezoidal rule would be preferable)
calc_area <- function(x, b, WSE)
{
  # N = discretization
  N <- length(x)
  delx <- (max(x)-min(x))/N
  A <- 0
  for (i in 1:N)
  {
    A <- A + (WSE-b[i])*delx
  }
  return(A)
}

wse.bf <- max(cross.sections[[r]]$b)
wse.bf <- b.min[1]+dbf[1] # not sure which it should be...
# A.riemann <- calc_area(cross.sections[[r]]$x, cross.sections[[r]]$b, wse.bf)



# --------------------------------------------------------------------------------------------------
# Wetted perimeter 

# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

library(pracma)

calc_area <- function(x, b, WSE)
{
  # Calculates area of irregular cross section using trapezoidal rule
  #
  # Inputs
  # x coordinates, bed elevation coordinates, WSE
  # Requires pracma package::trapz
  #
  # Outputs
  # Cross sectional area of flow
  
  keep.ind <- which(b<WSE)
  
  
  A.trap <- trapz(x, WSE-b)
  return(A.trap)
}

calc_WP <- function(x,b,WSE)
{
  # Calculates wetted perimeter of irregular cross section
  # For a shallow cross section, this is very close to its width
  #
  # Inputs
  # x coordinates, bed elevation coordinates, WSE
  #
  # Outputs
  # Wetted perimeter
  
  which(b<WSE)
  
  
  
  WP <- poly_length(x,b)
  
  return(WP)
}

# --------------------------------------------------------------------------------------------------

# function()





