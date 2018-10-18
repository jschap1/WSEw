#' Plot reach-averaged Pepsi Challenge data
#'
#' May need to edit some of the internal parameters here to get a good plot, e.g. ylim.
#' @examples plot_ra_pepsi_data(WSEw.unsorted, A)

plot_ra_pepsi_data <- function(WSEw.unsorted, A, reach_length = 10e3)
{

  par(mfrow = c(1,3))

  nr <- dim(A)[1]
  nt <- dim(A)[2]
  downstream_dist <- reach_length*(1:nr)
  t <- 1:nt

  # Area
  plot(t, A[1,], main = "Flow area",
       xlab = "time step (days)",
       ylab = "Flow area (sq. m)",
       ylim = c(0, max(A)),
       "n")
  cols <- rainbow(nr)
  for (r in 1:nr)
  {
    lines(t, A[r,], col = cols[r])
  }
  legend("bottomleft",
         legend = c("reach 1","reach 2","reach 3","reach 4"),
         col = cols, lty = c(1,1,1,1))

  # Height
  plot(t, WSEw.unsorted[[1]]$WSE, main = "Water surface elevation",
       xlab = "time step (days)",
       ylab = "h (m)",
       ylim = c(0, 20),
       "n")
  cols <- rainbow(nr)
  for (r in 1:nr)
  {
    lines(t, WSEw.unsorted[[r]]$WS, col = cols[r])
  }
  legend("bottomleft",
         legend = c("reach 1","reach 2","reach 3","reach 4"),
         col = cols, lty = c(1,1,1,1))

  # Width
  plot(t, WSEw.unsorted[[1]]$w, main = "Flow width",
       xlab = "time step (days)",
       ylab = "Flow width (m)",
       ylim = c(0, 1000),
       "n")
  cols <- rainbow(nr)
  for (r in 1:nr)
  {
    lines(t, WSEw.unsorted[[r]]$w, col = cols[r])
  }
  legend("bottomleft",
         legend = c("reach 1","reach 2","reach 3","reach 4"),
         col = cols, lty = c(1,1,1,1))

}
