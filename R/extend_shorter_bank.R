#' Make depths at left and right banks zero
#' 
#' Assumes the first and last point in the cross section are zero depth (banks)
#' @export
#' @param d list of depths at cross sections
#' @details May want to use a different assumption in the future

make_zero_bank <- function(d)
{
  nseg <- length(d)
  channel.pix <- unlist(lapply(d, length))
  for (seg in 1:nseg)
  {
    if (all(is.na(d[[seg]])))
    {
      next
    }
    d[[seg]][1] <- 0
    d[[seg]][channel.pix[seg]] <- 0
  }
  return(d)
}
