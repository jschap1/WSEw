#' Updates the WSEw package
#' 
#' Run from the /Users/.../Cross_Sections directory with devtools and roxygen2 packages loaded
#' @export
#' @details Sometimes you need to restart the R session in RStudio for the documentation to load properly
#' @example update_WSEw()

update_WSEw <- function()
{
  setwd("Codes/pkg")
  setwd("WSEw")
  document()
  setwd("..")
  install("WSEw")
  setwd("../..")
}

