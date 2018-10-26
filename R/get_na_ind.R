#' Get NA index
#' 
#' Gets the index of cross sections that are not to be included in the calculations
#' @export
# Some cross sections are just NA values.
# This is because the transect did not overlap any bathymetry data
# or because there are multiple channels.
get_na_ind <- function(cross_sections)
{
  n.xs <- length(cross_sections$x)
  na.ind <- vector(length = n.xs)
  for (i in 1:n.xs)
  {
    na.ind[i] <- all(is.na(cross_sections$d[[i]])) 
  }
  return(na.ind)
}