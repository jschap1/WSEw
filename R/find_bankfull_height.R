#' Find bankfull height
#' 
#' Finds the height at which overbank flow occurs
#' @param x distance
#' @param b bed elevation
#' @param delh increment for h
#' @details Takes (x,b) data at a cross section as input 
#' Works fairly well, but fails when a cross section has multiple subchannels
#' @export
#' @examples x = cross_sections$x[[7365]]
#' b = cross_sections$b[[7365]]
#' hbf <- find_bankfull_height(x,b)
#' plot(x, b, type = "l")
#' abline(hbf, 0)
#' @importFrom rootSolve uniroot.all

find_bankfull_height <- function(x, b, delh = 0.01)
{
  
  x.max <- max(x)
  
  # h <- min(b) # start at the minimum elevation
  h <- min(b) + 0.1*(max(b) - min(b)) # start at 10% of the bankfull depth to avoid issues with the jagged bottom
  
  f1 <- approxfun(x = x, y = b, method="linear", 
                  yleft = 1, yright = 1, rule = 1, f = 0, ties = mean)
  # as soon as there is only one intersection, then the bankfull height has been found
  first_iter <- TRUE
  xl_old <- -Inf
  xr_old <- Inf
  
  while (1)
  {

    h <- h + delh
    
    realistic_wd <- function(x, d = h)
    {
      pred <- f1(x)-d
      return(pred)
    }
    
    rts <- tryCatch(
      
      uniroot.all(realistic_wd, 
                  interval = c(0, x.max)
      ), 
      
      error = function(e) {NA}
      
    )
    
    if (length(rts) <= 1)
    {
      print("Only one intersection. Found hbf.")
      break
    }
    
    if (first_iter)
    {
      xl_old <- Inf
      xr_old <- -Inf
      first_iter <- FALSE
    } else
    {
      xl_old <- xl
      xr_old <- xr
    }
    
    # Find the left and right intersections
    gap <- which.max(diff(rts))
    xl <- rts[gap] 
    xr <- rts[gap+1]
    
    if ((xl > xl_old) | (xr < xr_old))
    {
      print("Found hbf.")
      break
    }
    
    # if (h>300) # to avoid infinite loop
    # {
    #   print("Reached maximum number of iterations. Did not find hbf.")
    #   break
    # }
     
  }
  
  bf <- list(hbf = h, xl = xl, xr = xr)
  
  return(bf)
  
}


# ----------------------------------------------------------------------------------------------------

# Minimum width-to-depth ratio method for identifying bankfull height
#
# cross_sections <- readRDS("./True_Parameters/p4/cross_sections.rds")
# xWSEw <- readRDS("./True_Parameters/p4/xWSEw.rds")
# 
# k <- 1709
# par(mfrow=c(2,2))
# plot(cross_sections$x[[k]], cross_sections$b[[k]], type = "l", col="red")
# plot(WSE~w, xWSEw[[k]])
# min_w_d(xWSEw[[k]])
# 
# cross_sections$n.channels # try a cross section with multiple channels
# 
# min_w_d <- function(WSEw)
# {
#   
#   d <- WSEw$WSE - min(WSEw$WSE)
#   d <- d[-1]
#   w <- WSEw$w
#   w <- w[-1]
#   plot(d + min(WSEw$WSE), w/d)
#   
#   which.min(w/d)
#   
#   
# }

  