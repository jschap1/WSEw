#' Plot time series
#' 
#' Plots a time series from a beginning date to an ending date
#' @param t Date object
#' @param x a variable 
#' @param begin.date date object telling the function on what date to start plotting the time series
#' @param end.date date object telling the function on what date to stop plotting the time series
#' @param ... other arguments as taken by plot()
#' @examples 
#' setwd("/Users/jschap/Documents/Research/SWOTBATH/")
#' load("./Data/USACE/Stage/stage_height_quincy.rda")
#' begin.date <- as.Date("2008-10-1") # beginning of water year
#' end.date <- as.Date("2011-9-30") # end of water year
#' my_plot_ts(t, h, begin.date, end.date, type = "l", 
#'            main = "Mississippi River at Quincy (high)", 
#'            xlab = "Date", ylab = "stage (m)",
#'            ylim = c(140, 149))
#' @details Really, it is preferable to use plot.ts()

my_plot_ts <- function(t, h, begin.date, end.date, ...)
{
  i <- which(t==begin.date)
  f <- which(t==end.date)
  
  plot(t[i:f], h[i:f], ...)
}



