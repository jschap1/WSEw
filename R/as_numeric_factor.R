#' As numeric factor
#'
#' Converts a factor to numeric, without losing information
#' @export

as_numeric_factor <- function(x) 
{
  as.numeric(levels(x))[x]
}
