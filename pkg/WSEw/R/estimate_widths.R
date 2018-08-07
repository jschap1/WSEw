#' Estimate Widths
#' 
#' Estimates bankfull widths of transects based on the angle of the transect and the number of pixels in the transect.
#' A simplification. Hopefully, it is accurate over large enough transects. 
#' Should probably revisit and make sure this code is reasonably accurate.
#' @param rpolyline a polyline crossing the river channel
#' @param tlength 
#' @param nseg 
#' @param resolution resolution of bathymetry data (m)
#' @keywords bankfull width
#' @import sp
#' @export
#' @examples 
#' wbf <- estimate_widths(rpolyline, resolution = 5, channel.pix, nseg)

estimate_widths <- function(rpolyline, tlength, nseg, resolution = 5)
{
  cross_section_width <- vector(length=nseg)
  for (seg in 1:nseg)
  {
    theta <- atan((rpolyline$y[seg+1] - rpolyline$y[seg])/(rpolyline$x[seg+1] - rpolyline$x[seg]))
    if(is.na(theta)) {theta <- pi/2} # in case of vertical segment
    theta_perp <- theta + pi/2
    cross_section_width[seg] <- tlength[[seg]]*abs(cos(theta_perp))*resolution
  }
  return(cross_section_width)
}