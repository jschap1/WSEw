#' Hydraulic parameter predictions boxplots
#' 
#' Makes box and whisker plots showing predicted hydraulic parameter error
#' @param z.l linear model parameter predictions or prediction errors
#' @param z.sb SB model parameter predictions or prediction errors
#' @param z.nl NL model parameter predictions or prediction errors
#' @param z.nlsb NLSB model parameter predictions or prediction errors
#' @param k exposure level
#' @param lumped flag, if true, lump the cross sections together and just compare methods
#' @param legend flag for legend
#' @param absolute flag for plotting absolute error
#' @param ... other parameters as taken by boxplot()
#' @export
#' @details If lumped = FALSE, plots the results over each cross section separately.
#' This is the updated version, which can take either predictions or prediction error as input.
#' @examples 
#' Test case
#' z.l <- z0.l.error[1:3,,1:2]
#' z.sb <- z0.sb.error[1:3,,1:2]
#' z.nl <- z0.nl.error[1:3,,1:2]
#' z.nlsb <- z0.nlsb.error[1:3,,1:2]
#' parameter_predictions_boxplots(z.l, z.sb, z.nl, z.nlsb, k = 12, names = bplab)
#' z.l[,12,]
#' z.sb[,12,]
#' z.nl[,12,]
#' z.nlsb[,12,]
#' 
#' bplab <- c("1","1","1","1","2","2","2","2","3","3","3","3")
#' lumpedlab <- c("L","SB","NL","NLSB")
#' par(mfrow = c(2,2))
#' par(opar)
#' k <- 16
#' z0 error
#' parameter_predictions_boxplots(z0.l.error, z0.sb.error, z0.nl.error, z0.nlsb.error,
#'                                k = k,
#'                                lumped = TRUE,
#'                                main = paste("z0 error at", 100*expo[k], "percent channel exposure"),
#'                                ylab = "z0 error (m)",
#'                                legend = FALSE,
#'                                notch = TRUE,
#'                                names = lumpedlab)

#' A0 error
#' parameter_predictions_boxplots(A0.l.error, A0.sb.error, A0.nl.error, A0.nlsb.error,
#'                                k = k,
#'                                lumped = TRUE,
#'                                main = paste("A0 error at", 100*expo[k], "percent channel exposure"),
#'                                ylab = "A0 error (sq. m)",
#'                                legend = FALSE,
#'                                notch = TRUE,
#'                                names = lumpedlab)
#' 
#' parameter_predictions_boxplots(z0.l, z0.sb, z0.nl, z0.nlsb, z0.true.ra, k = k, lumped = FALSE,
#'                                main = paste("z0 error at", 100*expo[k], "percent channel exposure"),
#'                                ylab = "average z0 error (m)", legend = FALSE, A0 = FALSE, notch = TRUE, names = bplab)
#' abline(h = z0.true.ra[3])

parameter_predictions_boxplots <- function(z.l, z.sb, z.nl, z.nlsb, k,
                                           lumped = FALSE, legend = TRUE, ...)
{
 
  nr <-  dim(z.l)[1] # get number of reach-average cross sections
  M <- dim(z.l)[3] # get number of replicates
  
  z <- c(as.vector(t(z.l[,k,])), # (as.vector makes a matrix into a vector by column)
         as.vector(t(z.sb[,k,])),
         as.vector(t(z.nl[,k,])),
         as.vector(t(z.nlsb[,k,])))
  
  xs <- rep(rep(1:nr, each = M), 4)
  
  type <- c(rep("L", nr*M), rep("SB", nr*M), rep("NL", nr*M), rep("NLSB", nr*M))
  
  df <- data.frame(z = z, xs = xs, type = type)
 
  df$type=factor(df$type, levels=levels(df$type)[c(1,4,2,3)]) # reorder
  
  cols <- c("red","orange","green","blue")
  
  if (lumped)
  {
    b <- boxplot(z ~ type, df, col = cols, ...)
    abline(0,0, lty = 2, lwd = 1.5)
    
  } else
  {
    b <- boxplot(z ~ type + xs, df, col = cols, ...)
    abline(0,0, lty = 2, lwd = 1.5)
    if (legend)
    {
      legend("topleft", horiz = TRUE, legend = c("L","SB","NL","NLSB"), fill = cols)
    }
  }
  
  return(b)
}
