#' Slope Break Fit for WSE-w Relationship
#' 
#' Optionally uses Mersel's original method, which screens out "non-optimal" cross sections
#' @param WSEw WSEw data (at a given level of exposure)
#' @param mersel flag for Mersel method
#' @param thres threshold for Mersel method, must be provided if using Mersel method. 0.015 has been used in the past.
#' @param multiple_breaks one slope break or multiple slope breaks
#' @param continuity require continuity between piecewise fits?
#' @export

fit_slopebreak <- function(WSEw,  mersel = FALSE, thres = NULL, multiple_breaks = FALSE, continuity = TRUE)
{
  
  nn <- length(WSEw$w) # number of data points
  
  if (mersel)
  {
    # Use original Mersel method
    
    
    
    
  } else
  {
    
    if (multiple_breaks)
    {
      # Use multiple slope breaks
      b <- breakpoints(WSE~w, data = WSEw, h=5)$breakpoints # h is the minimum number of points required for a section
      if (is.null(b)) {b<-NA}
      
      
      
    } else
    {
      # Use one slope break
      b <- breakpoints(WSE~w, data = WSEw, breaks = 1, h=5)$breakpoints
      if (is.null(b)) 
      {
        b<-NA
        print("No breakpoints could be identified")
      }
      
      sb.ind <- b[1]
      WSEw1 <- WSEw[1:sb.ind,]
      WSEw2 <- WSEw[sb.ind:nn,]
      
      if (!continuity)
      {
        lf1 <- lm(WSE~w, data = WSEw1) # below slope break
        lf2 <- lm(WSE~w, data = WSEw2) # above slope break
      } else if (continuity)
      {
        lf1 <- lm(WSE~w, data = WSEw1) # might be better off replacing this with a built-in function, like the segmented pkg
        a0 <- as.numeric(coef(lf1)[1])
        a1 <- as.numeric(coef(lf1)[2])
        wb <- WSEw$w[sb.ind]
        intercept <- a0+a1*wb
        lf2 <- lm(WSE ~ -1 + I(w-wb), data = WSEw2, offset = rep(intercept,dim(WSEw2)[1]))
      }
      fits <- list(lf1, lf2)
      attributes(fits) <- list(sb.ind = sb.ind) # adding the sb.ind as an output
    }
  }
  
  return(fits)
}


