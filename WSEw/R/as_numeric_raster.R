#' As Numeric Raster
#' 
#' Converts a raster of factor values to the a raster of the corresponding numeric values.
#' Somewhat more efficient than get_depth_from_lutable, but can still take a long time for big rasters.
#' See raster::factor documentation.
#' @export
#' @param r a raster with an attribute table
#' @param att column of the attribute to extract
#' @details Given a raster with data values coded in an attribute table, creates a raster directly composed of those values. 
#' Loops over every entry in the raster; it is not very efficient. 
#' @return r2 raster of numeric values coded from the input raster
#' @examples
#' umesc <- raster("/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC/bath_pool_21/bath_1999_p21")
#' umesc <- raster("bath_pool_21/bath_1999_p21")
#' depth <- as_numeric_raster(umesc, att = "DEPTH_M")
#' @keywords deratify factor raster
#' @importFrom raster deratify

as_numeric_raster <- function(r, att)
{
  
  lutable <- levels(r)[[1]]
  r1 <- r
  # r1 <- deratify(r, att = att) # this step does not appear to be necessary
  data.ind <- which(!is.na(getValues(r1))) # only operates over cells with data to save time
  factor.vals <- r1[data.ind]
  
  ##################################################
  N <- length(factor.vals)
  numeric.vals <- vector(length = N)
  for (k in 1:N)
  {
    # currently only set up for "DEPTH_M" attribute
    numeric.vals[k] <- as.numeric.factor(lutable$DEPTH_M[match(factor.vals[k],lutable$ID)])
    # Problem: factor.vals is not necessarily in lutable$ID
    if (k%%50000 == 0)
    { 
      print(paste0("progress: ", round(100*k/N, 2), "%")) # display progress
    }
  }
  ##################################################

  r2 <- r1
  r2[data.ind] <- numeric.vals
  plot(r2)

  return(r2)
}

# Helper function
as.numeric.factor <- function(x) 
{
  as.numeric(levels(x))[x]
}
