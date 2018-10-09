#' Get Depth from Lookup Table
#' 
#' Not used as of 8/3/2018 update. Replaced by as_numeric_raster.
#' Extracts depth from a lookup table stored with the bathymetry raster. 
#' The UMESC bathymetry data are stored using a special format. 
#' This function decodes that special format and puts the data in units of meters.
#' @param transect depth values along a transect
#' @keywords lookup table, data compression
#' @export
#' @examples 
#' get_depth_from_lutable(transect)
 
get_depth_from_lutable <- function(transect, depth)
{
  lutable <- depth@data@attributes[[1]]
  transect.depth <- vector(length = length(transect))
  for (k in 1:length(transect))
  {
    transect.depth[k] <- as.numeric.factor(lutable$DEPTH_M[match(transect[k],lutable$ID)])
    if (k%%50000 == 0)
    { 
      print(k) # display progress
    }
  }
  return(transect.depth)
}

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}




