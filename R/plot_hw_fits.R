#' Plot h-w fits for one cross section
#' 
#' Use exposure = 1 for Pepsi Challenge data. Can optionally add measurement error.
#' @param ... other parameters as taken by plot()
#' @return predicted hydraulic parameters
#' @export
#' @examples r <- 1
#' pred <- plot_hw_fits(rWSEw[[r]], exposure = 1, 
#'                     xlim = c(0, 600), ylim = c(131,148), 
#'                     main = paste("Reach", r),
#'                     xlab = "w (m)", ylab = "WSE (m)")

plot_hw_fits <- function(WSEw, exposure, sd_wse = 0, sd_w = 0, ...)
{
  
  WSEw_obs <- observe(WSEw, exposure = exposure, sd_wse = sd_wse, sd_w = sd_w)
  
  # Fit models
  lf <- fit_linear(WSEw_obs)
  sb <- fit_slopebreak(WSEw_obs, multiple_breaks = FALSE, continuity = TRUE)
  nl <- fit_nonlinear(WSEw_obs)
  nlsb <- fit_nlsb(WSEw_obs)
  
  plot(WSE~w, WSEw_obs, ...) # observations, with error
  lines(WSE~w, WSEw) # true data, without error
  
  # plot linear
  lines(WSEw_obs$w, predict(lf), col = "red")
  w.vals <- seq(0, min(WSEw_obs$w), 1)
  lines(w.vals, predict(lf, newdata = data.frame(w=w.vals)), 
        col = "red", lty = 2)
  
  # plot slope break
  nn <- length(WSEw_obs$w)
  sb.ind <- attributes(sb)$sb.ind
  lines(WSEw_obs$w[1:sb.ind], predict(sb[[1]]), col = "orange")
  lines(WSEw_obs$w[(sb.ind):nn], predict(sb[[2]]), col = "orange")
  lines(w.vals, predict(sb[[1]], newdata = data.frame(w=w.vals)), 
        col = "orange", lty = 2)
  
  # plot nl
  lines(WSEw_obs$w, predict(nl), col = "green")
  lines(w.vals, predict(nl, newdata = data.frame(w=w.vals)), 
        col = "green", lty = 2)
  
  # plot nlsb
  nlsb.ind <- attributes(nlsb)$sb.ind
  lines(WSEw_obs$w[1:nlsb.ind], predict(nlsb[[1]]), col = "blue")
  lines(WSEw_obs$w[(nlsb.ind):nn], predict(nlsb[[2]]), col = "blue")
  lines(w.vals, predict(nlsb[[1]], newdata = data.frame(w=w.vals)), 
        col = "blue", lty = 2)

  # z0 predictions
  z0.l <- as.numeric(coef(lf)[1])
  z0.sb <- as.numeric(coef(sb[[1]])[1])
  z0.nl <- as.numeric(coef(nl)[1])
  z0.nlsb <- as.numeric(coef(nlsb[[1]])[1])
  
  # A0 predictions
  A0.l <- calc_model_A0(lf, type = "linear", pos.only = FALSE)
  A0.sb <- calc_model_A0(sb, type = "sb", pos.only = FALSE)
  w1 <- min(WSEw_obs$w)
  h1 <- min(WSEw_obs$WSE)
  A0.nl <- calc_model_A0(nl, type = "nl", w1 = w1, h1 = h1, pos.only = FALSE)
  A0.nlsb <- calc_model_A0(nlsb, type = "nlsb", w1 = w1, h1 = h1, pos.only = FALSE)
  
  pred <- vector(length = 2, "list")
  pred$z0 <- list(z0.l = z0.l, z0.sb = z0.sb, z0.nl = z0.nl, z0.nlsb = z0.nlsb)
  pred$A0 <- list(A0.l = A0.l, A0.sb = A0.sb, A0.nl = A0.nl, A0.nlsb = A0.nlsb)
  return(pred)
}
