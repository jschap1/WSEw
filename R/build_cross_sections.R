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
  dbf <- unlist(lapply(d, max))
  delb <- vector(length = nseg) # initialize, only used if option = zero, though
  
  if (option == "zero")
  {
    
    # save the original depths
    dl <- vector(length = nseg)
    dr <- vector(length = nseg)
    for (seg in 1:nseg)
    {
      dl[seg] <- d[[seg]][1]
      dr[seg] <- d[[seg]][length(d[[seg]])]
      
      # how large the adjustment to 0 is, as a function of bankfull depth
      delb[seg] <- max(dl[seg], dr[seg])/dbf[seg] 
    }
    
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
  
  if (option == "truncate")
  {
    
    for (seg in 1:nseg)
    {
      if (channel.pix[seg] <= 1)
      {
        next
      }
      bf <- find_bankfull_height(x[[seg]],b[[seg]])
      l.ind <- which.min(abs(x[[seg]]-bf$xl))
      r.ind <- which.min(abs(x[[seg]]-bf$xr))
      b[[seg]] <- b[[seg]][l.ind:r.ind]
      x[[seg]] <- seq(0, x[[seg]][r.ind] - x[[seg]][l.ind], length.out = length(b[[seg]]))
      d[[seg]] <- refWSE - b[[seg]]
    }

  }

  cross.sections <- list(x = x, b = b, d = d, channel.pix = channel.pix, delb = delb)
  return(cross.sections)
}