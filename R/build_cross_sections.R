#' Build cross sections
#'
#' Stores cross section data (bed elevation, depth, and x coordinate) in a convenient format
#' @export
#' @param option Use "zero" to make the left and right banks have depths zero
#' Use "extend" to extend the shorter bank height to match the taller bank height
#' Use "truncate" to truncate the taller bank height to match the shorter bank height.

build_cross_sections <- function(d, refWSE, wbf, option = "none")
{
  # Get number of values in each transect
  # This is done to screen out transects that are unrealistically narrow
  channel.pix <- unlist(lapply(d, length))
  nseg <- length(d)

  if (option == "extend")
  {
    d <- extend_shorter_bank(d)
  } else if (option == "truncate")
  {
    d <- truncate_taller_bank(d)
  } else if (option == "zero")
  {
    d <- make_zero_bank(d)
  }
  
  x <- vector("list", length = nseg)
  for (seg in 1:nseg)
  {
    if (all(is.na(d[[seg]])))
    {
      next
    }
    x[[seg]] <- seq(0, wbf[seg], length.out = length(d[[seg]]))
  }

  b <- lapply(d, function(x, refWSE) {refWSE-x}, refWSE)
  cross.sections <- list(x = x, b = b, d = d, channel.pix = channel.pix)
  return(cross.sections)
}