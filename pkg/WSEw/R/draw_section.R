#' Draw cross section
#'
#' Draws the modeled cross section superimposed on the true cross section
#' Cannot actually use "reach-averaged cross sections" because it is unclear how to define them
#' @param r reach number
#' @param x cross section number, should be roughly representative of the reach-avg
#' @param cross_sections cross section geometry
#' @param model model
#' @param type
#' @export
#' 
#' In development.

draw_section <- function(r, x, cross_sections, model, type)
{
  
  plot(cross_sections$x[[x]], cross_sections$b[[x]], type = "l", 
       xlab = "x (m)", ylab = "b (m)")
  
  if (type == "linear")
  {
    
    
    
  } else if (type  "sb")
  {
    
  } else if (type  "sbm")
  {
    
  } else if (type  "nl")
  {
    
  } else if (type  "nlsb")
  {
    
  }
  
  
  
  
}