#' Count separate channels
#' 
#' count number of separate channels
#' @param t = transect
#' @export

count_separate_channels <- function(t)
{
  changeflag <- FALSE
  prev.val <- NA
  n <- length(t)
  num.channels <- 0
  for (i in 1:n)
  {
    if (!is.na(t[i]) && is.na(prev.val))
    {
      changeflag <- TRUE
    }
    if (changeflag)
    {
      num.channels <- num.channels + 1
    } 
    prev.val <- t[i]
    changeflag <- FALSE
  }
  return(num.channels)
}

# ------------------------------------------------------------------------------------------------
#' Count separate channels 2
#' 
#' count number of separate channels, can filter by min_width
#' @details Can adjust min_width to exclude very narrow channels
#' @param t = transect
#' @export

count_separate_channels_2 <- function(t, min_width = 1)
{
  w <- 0
  changeflag <- FALSE
  prev.val <- NA
  n <- length(t)
  num.channels <- 0
  for (i in 1:n)
  {
    if (is.na(t[i]) && is.na(prev.val))
    {
      w <- 0
    } else if (!is.na(t[i]) && is.na(prev.val))
    {
      w <- 1
      if (i == n && w >= min_width)
      {
        num.channels <- num.channels + 1
      }
    } else if (!is.na(t[i]) && !is.na(prev.val))
    {
      w <- w + 1
    } else if (is.na(t[i]) && !is.na(prev.val))
    {
      # print(w)
      if (w >= min_width)
      {
        num.channels <- num.channels + 1
      }
      w <- 0
    }
    prev.val <- t[i]
    changeflag <- FALSE
  }
  return(num.channels)
}


  



