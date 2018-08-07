#' As Numeric Raster
#' 
#' Converts a raster of factor values to the a raster of the corresponding numeric values.
#' Somewhat more efficient than get_depth_from_lutable, but can still take a long time for big rasters.
#' See raster::factor documentation
#' @param r a raster with an attribute table
#' @param att column of the attribute to extract
#' @keywords deratify factor raster
#' @importFrom raster deratify
#' @export
as_numeric_raster <- function(r, att)
{
  
  lutable <- levels(r)[[1]]
  r1 <- deratify(r, att = att) # this step might not be necessary
  data.ind <- which(!is.na(getValues(r1))) # only operates over cells with data to save time
  factor.vals <- r1[data.ind]
  
  ##################################################
  N <- length(factor.vals)
  numeric.vals <- vector(length = N)
  for (k in 1:N)
  {
    # currently only set up for "DEPTH_M" attribute
    numeric.vals[k] <- as.numeric.factor(lutable$DEPTH_M[match(factor.vals[k],lutable$ID)]) 
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
