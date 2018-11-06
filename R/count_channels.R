#' Count channels
#' 
#' Returns the number of channels contained within a (possibly multi-channel) cross section
#' @param x x coordinate
#' @param b bed elevation
#' @details Finds the local minima in the channel cross section
#' Uses smoothing techniques described at https://rpubs.com/mengxu/peak_detection
#' Accuracy depends on choice of tuning parameters
#' Could calibrate to a particular dataset
#' @examples 
#' k <- 1000
#' x <- cross_sections$x[[k]]
#' b <- cross_sections$b[[k]]
#' d <- cross_sections$d[[k]]
#' plot(x,b,type="l")
#' @export

count_channels <- function(x, d)
{
  
  # skip empty cross sections
  peaks <- test(x = x, y = d, w = 5, span = 0.55, plotflag = FALSE)
  nchannels <- length(peaks$x)
  return(nchannels)
  
}

# ----------------------------------------------------------------------------------------------------

#' Test
#' 
#' Called Test because it's supposed to be a way to find good tuning parameters using trial and error
#' @param w tuning parameter
#' @param span tuning parameter

test <- function(x, y, w, span, plotflag = TRUE) 
{
  
  peaks <- argmax(x, y, w=w, span=span)
  
  if(plotflag)
  {
    plot(x, y, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep=""))
    lines(x, peaks$y.hat,  lwd=2) #$
    y.min <- min(y)
    sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, peaks$y.hat[i]), col="Red", lty=2))
    points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, cex=1.25)
  }

  return(peaks)
  
}

# ----------------------------------------------------------------------------------------------------

argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}
