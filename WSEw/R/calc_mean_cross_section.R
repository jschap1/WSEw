#' Calculate Mean Cross Section
#' 
#' Uses linear interpolation to calculate the average cross section from 
#' among individual cross sections with different lengths
#' @param xs.res cross section geometry that you wish to average, 
#' resampled to the same length 
#' @export
#' @examples xs.avg <- calc_mean_cross_section(xs.res, n = 100)
#' calc_mean_cross_section(xs.res[1:100,,])

calc_mean_cross_section <- function(xs.res)
{
  mean.x <- apply(xs.res[,,1], 2, mean)
  mean.b <- apply(xs.res[,,2], 2, mean)
  mean.d <- apply(xs.res[,,3], 2, mean)
  xs.avg <- data.frame(x = mean.x, b = mean.b, d = mean.d)
  return(xs.avg)
}

#' Resample Cross Section
#' 
#' Uses linear interpolation to resample cross sections with different resolutions
#'to the same resolution.
#' @param cross_sections list of cross section data  
#' @param n resolution of the resample cross sections
#' @examples
#' xs.res <- resample_xs(cross_sections, n)
#' # check that it worked
#' x <- 100
#' plot(cross_sections$x[[x]], cross_sections$b[[x]])
#' lines(xs.res[x,,1], xs.res[x,,2], col = "red")
#' @export

resample_xs <- function(cross_sections, n)
{
  n.xs <- length(cross_sections$x)
  
  # new data structure where all cross sections have the same number of points
  # data frame would be a better output format than array bc names
  xs.res <- array(dim = c(n.xs, n, 3)) # cross sections by numdata by (x,b,d)
  
  for(x in 1:n.xs) # for each cross section
  {
    xs <- data.frame(x = cross_sections$x[[x]], b = cross_sections$b[[x]]) # get coordinates of cross section
    d.coord <- cross_sections$d[[x]]
    xs.res[x,,] <- cbind(x = approx(xs, n = n)$x,  # linearly interpolate to the new length n
                         b = approx(xs, n = n)$y, 
                         d = approx(d.coord, n = n)$y)
  }
  return(xs.res)
  
}







