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
#' 
#' Here is some code I've been playing with:
#' pool <- 21
#' i <- 8
#' d <- readRDS(paste0("./Outputs/p", pool, "/main_channel.rds"))
#' wbf <- readRDS(paste0("./Outputs/p", pool, "/wbf.rds"))
#' cross_sections <- build_cross_sections(d, refWSE[i], wbf, option = "zero")
#' head(which(cross_sections$channel.pix>100))
#' na.ind <- get_na_ind(cross_sections)
#' head(which(!na.ind), 100)
#' x <- 4041
#' plot(cross_sections$x[[x]], cross_sections$b[[x]], type = "l")
#' count_channels(cross_sections$x[[x]], cross_sections$d[[x]])
#' 
#' w1 <- round(0.2*length(cross_sections$x[[x]])) #' one idea for the w parameter's value
#' 100/wbf[x]
#' 
#' #' Make an animation of how the river channel changes from upstream to downstream
#' na.ind <- get_na_ind(cross_sections)
#' for (x in which(!na.ind))
#' {
#'   num <- sprintf("%05d", x) #' ensuring the same number of digits so they files can be sorted in the proper order
#'   png(paste0("pool21_xs_", num, ".png"))
#'   plot(cross_sections$x[[x]], 
#'        cross_sections$b[[x]], 
#'        type = "l", 
#'        main = paste("xs", x, "wbf = ", round(wbf[x]))
#'   )
#'   dev.off()
#'   #' readline(prompt="Press [enter] to continue")
#' }
#' 
#' nc <- vector(length = length(which(!na.ind)))
#' for (x in which(!na.ind))
#' {
#'   nc[x] <- count_channels(cross_sections$x[[x]], cross_sections$d[[x]])
#'   if (x>2000){break}
#' }


#' 
#' @export

count_channels <- function(x, d)
{
  wbf <- max(x)
  swot.res <- 100
  # w1 <- min(round(0.1*length(x)), 20)
  w1 <- round(0.1*length(x))
  peaks <- test(x = x, y = d, w = w1, span = swot.res/wbf, plotflag = FALSE)
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
