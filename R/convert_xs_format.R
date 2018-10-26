#' Convert cross section format
#'
#' Converts cross section data structure to a list of data frames'
convert_xs_format <- function(x, b, d)
{
  nseg <- length(x)
  cross.sections <- vector(length = nseg, "list")
  for (seg in 1:nseg)
  {
    cross.sections[[seg]] <- data.frame(x = x[[seg]], b = b[[seg]], d = d[[seg]])
  }
  return(cross.sections)
}
