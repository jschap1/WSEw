#' Updates the WSEw package
#' 
#' Run from the /Users/.../Cross_Sections directory with devtools and roxygen2 packages loaded
#' @export

update_WSEw <- function()
{
  setwd("Codes/pkg")
  setwd("WSEw")
  document()
  setwd("..")
  install("WSEw")
  setwd("../..")
}

