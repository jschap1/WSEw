#' Build cross sections
#'
#' Stores cross section data (bed elevation, depth, and x coordinate) in a convenient format
#' @export

build_cross_sections <- function(d, refWSE, wbf)
{
  # Get number of values in each transect
  # This is done to screen out transects that are unrealistically narrow
  channel.pix <- unlist(lapply(d, length))

  # Assume the first and last point are zero depth (banks)
  # May want to use a different assumption in the future
  nseg <- length(wbf)
  for (seg in 1:nseg)
  {
    if (all(is.na(d[[seg]])))
    {
      next
    }
    d[[seg]][1] <- 0
    d[[seg]][channel.pix[seg]] <- 0
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

  # Remove (nearly) empty channels
  # This may not be necessary
  # if (any(channel.pix<=10))
  # {
  #   rm.ind <- which(channel.pix<=10)
  #   d <- d[-rm.ind]
  #   b <- b[-rm.ind]
  #   x <- x[-rm.ind]
  # }

  cross.sections <- list(x = x, b = b, d = d, channel.pix = channel.pix)
  return(cross.sections)
}
