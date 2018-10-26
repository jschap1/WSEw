#' Smooth cross sections
#'
#' Smooths bed elevation and depth data using a moving average
smooth_cross_sections <- function(cross_sections, k)
{

  nseg <- length(cross_sections$x)

  b.smooth <- vector("list", length = nseg)
  d.smooth <- vector("list", length = nseg)
  for (seg in 1:nseg)
  {
    if (channel.pix[seg] <= k)
    { # error handling for zero-width channels
      next
    }
    b.smooth[[seg]] <- filter(b[[seg]], sides = 2, filter = rep(1/k,k))
    d.smooth[[seg]] <- filter(d[[seg]], sides = 2, filter = rep(1/k,k))
  }

  cross_sections <- cbind(cross_sections,
    b.smooth = b.smooth,
    d.smooth = d.smooth)

  return(cross_sections)

}



# ----------------------------------------------------------------------------------------------------
# Scrap

#
# # Remove null (empty) transects
# null.ind <- unlist(lapply(main_channel, is.null))
# na.ind <- lapply(main_channel, is.na)
# na.ind <- unlist(lapply(na.ind, all))
# if (sum(null.ind)>0)
# {
#   main_channel <- main_channel[-which(null.ind)]
#   xs_locs <- xs_locs[-which(null.ind)]
# }
# if (sum(na.ind)>0)
# {
#   main_channel <- main_channel[-which(na.ind)]
# }
# nseg <- length(main_channel)
#
# if(any(is.na(wbf)))
# {
#   wbf.na.ind <- which(is.na(wbf))
#   wbf <- wbf[-wbf.na.ind]
# }

# plot(main_channel[[1]], type = "l")
# I suspect the scaling factor is 1000, but double check.
