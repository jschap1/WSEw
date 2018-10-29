# Functions for modifying the bank heights (for computational reasons)
# This may not be the right way to do this because the left and right banks 
# are not necessarily the tallest
# Try doing something like this later.

# ----------------------------------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------------------------------

#' Extend shorter bank
#' 
#' Assumes the first and last point in the cross section are zero depth (banks)
#' @export
#' @param d list of depths at cross sections
#' @details May want to use a different assumption in the future

extend_shorter_bank <- function(d)
{
  nseg <- length(d)
  channel.pix <- unlist(lapply(d, length))
  for (seg in 1:nseg)
  {
    if (all(is.na(d[[seg]])))
    {
      next
    }
    taller.ind <- which.max(c(d[[seg]][1], d[[seg]][channel.pix[seg]]))
    if (taller.ind == 2)
    {
      d[[seg]][channel.pix[seg]] <- d[[seg]][1]
    } else
    {
      d[[seg]][1] <- d[[seg]][channel.pix[seg]]
    }
    
  }
  return(d)
}


# ----------------------------------------------------------------------------------------------------

#' Truncate taller bank
#' 
#' Assumes the first and last point in the cross section are zero depth (banks)
#' @export
#' @param d list of depths at cross sections
#' @details May want to use a different assumption in the future

truncate_taller_bank <- function(d)
{
  nseg <- length(d)
  channel.pix <- unlist(lapply(d, length))
  for (seg in 1:nseg)
  {
    if (all(is.na(d[[seg]])))
    {
      next
    }
    taller.ind <- which.max(c(d[[seg]][1], d[[seg]][channel.pix[seg]]))
    if (taller.ind == 2)
    {
      d[[seg]][1] <- d[[seg]][channel.pix[seg]]
    } else
    {
      d[[seg]][channel.pix[seg]] <- d[[seg]][1]
    }
    
  }
  return(d)
}