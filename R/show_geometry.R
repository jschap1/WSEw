#' Show geometry
#' 
#' Shows the cross section geometry and fitted model for a given reach and exposure level
#' @param WSEw true height and width values
#' @param WSEw_obs observed height and width values
#' @param cross_section data frame with x, b for a given reach
#' @param method l, sb, nl, nlsb
#' @param ... optional parameters for plot()
#' @examples
#' library(minpack.lm) # required for NL, NLSB methods
#' exp_dir <- "./Outputs/Final/p21" # pool 21
#' r <- 1 # reach 1
#' cross_sections_avg <- readRDS(file.path(exp_dir, "cross_sections_avg.rds"))
#' cross_section <- data.frame(x = cross_sections_avg$x[[r]], b = cross_sections_avg$b[[r]] )
#' rWSEw <- readRDS(file.path(exp_dir, "rWSEw.rds"))
#' WSEw_obs1 <- observe2(WSEw = rWSEw[[r]], exposure = 0.5)
#' show_geometry(WSEw = rWSEw[[r]], 
#'               WSEw_obs = WSEw_obs1, 
#'               cross_section = cross_section,
#'               method = "nlsb", ylim = c(130, 145))
#' @details The code could use a bit of tweaking, but it more or less does what it's supposed to.
#' @export
show_geometry <- function(WSEw, WSEw_obs, cross_section, method, showLegend = TRUE, ...)
{

  plot(cross_section$x, cross_section$b, 
       type = "l",
       xlab = "x (m)",
       ylab = "bed elevation (m)",
       col = "gray",
       ...)
  
  # ------------------------------------------------------------  
  # Linear
  
  if (method == "l")
  {
    lf1 <- fit_linear(WSEw_obs, h = 10)
    
    wbf <- max(cross_section$x)
    wmin <- min(WSEw_obs$w)
    y <- (wbf - wmin)/2
    
    x3 <- seq(wbf/2, wbf - y)
    x4 <- seq(wbf - y, wbf)
    
    b3 <- predict(lf1, newdata = data.frame(w = 2*(x3-wbf/2)))
    b4 <- predict(lf1, newdata = data.frame(w = 2*(x4-wbf/2)))
    
    lines(x3, b3, col = "red", lty = 2, lwd = 3)
    lines(x4, b4, col = "red", lty = 1, lwd = 3)
    
    b2 <- rev(b3)
    b1 <- rev(b4)
    
    x1 <- seq(0, y)
    x2 <- seq(y, wbf/2)
    lines(x1, b1, col = "red", lty = 1, lwd = 3)
    lines(x2, b2, col = "red", lty = 2, lwd = 3)
    
    if (showLegend)
    {
      legend("top", legend = c("True","Fitted","Extrapolated"), 
             col = c("gray", "red", "red"),
             lty = c(1,1,2))
    }

  }

  # ------------------------------------------------------------
  # SB
  
  if (method == "sb")
  {
    
    sb1 <- fit_slopebreak(WSEw_obs, multiple_breaks = FALSE, continuity = TRUE)
    
    wbf <- max(cross_section$x)
    wmin <- min(WSEw_obs$w)
    sb.ind <- attributes(sb1)$sb.ind
    wsb <- WSEw_obs$w[sb.ind]
    y <- (wbf - wsb)/2
    
    x4 <- seq(wbf/2, (wbf + wmin)/2)
    x5 <- seq((wbf + wmin)/2, wbf-y)
    x6 <- seq(wbf-y, wbf)
    
    b4 <- predict(sb1[[1]], newdata = data.frame(w = 2*(x4-wbf/2)))
    b5 <- predict(sb1[[1]], newdata = data.frame(w = 2*(x5-wbf/2)))
    # b6 <- predict(sb1[[2]], newdata = data.frame(w = 2*(x6-wbf/2)))
    
    slope <- as.numeric(coef(sb1[[2]])) # workaround for b6
    intercept <- predict(sb1[[1]], newdata = data.frame(w = wsb))
    predict_sb2 <- function(slope, intercept, x, offset = wsb)
    {
      pred <- intercept + slope*(x - offset)
    }
    b6 <- predict_sb2(slope, intercept, x = 2*(x6-wbf/2))
    
    x1 <- seq(0, y)
    x2 <- seq(y, (wbf - wmin)/2)
    x3 <- seq((wbf - wmin)/2, wbf/2)
    
    b1 <- rev(b6)
    b2 <- rev(b5)
    b3 <- rev(b4)

    lines(x1, b1, col = "orange", lty = 1, lwd = 3)
    lines(x2, b2, col = "orange", lty = 1, lwd = 3)
    lines(x3, b3, col = "orange", lty = 2, lwd = 3)
    lines(x4, b4, col = "orange", lty = 2, lwd = 3)
    lines(x5, b5, col = "orange", lty = 1, lwd = 3)
    lines(x6, b6, col = "orange", lty = 1, lwd = 3)
        
    if (showLegend)
    {
      legend("top", legend = c("True","Fitted","Extrapolated"), 
             col = c("gray", "orange", "orange"),
             lty = c(1,1,2))
    }
    
  }

  # ------------------------------------------------------------  
  # NL
  
  if (method == "nl")
  {
    nl1 <- fit_nonlinear(WSEw_obs)
    
    wbf <- max(cross_section$x)
    wmin <- min(WSEw_obs$w)
    y <- (wbf - wmin)/2
    
    x3 <- seq(wbf/2, wbf - y)
    x4 <- seq(wbf - y, wbf)
    
    b3 <- predict(nl1, newdata = data.frame(w = 2*(x3-wbf/2)))
    b4 <- predict(nl1, newdata = data.frame(w = 2*(x4-wbf/2)))
    
    lines(x3, b3, col = "green", lty = 2, lwd = 3)
    lines(x4, b4, col = "green", lty = 1, lwd = 3)
    
    b2 <- rev(b3)
    b1 <- rev(b4)
    
    x1 <- seq(0, y)
    x2 <- seq(y, wbf/2)
    lines(x1, b1, col = "green", lty = 1, lwd = 3)
    lines(x2, b2, col = "green", lty = 2, lwd = 3)
    
    if (showLegend)
    {
      legend("top", legend = c("True","Fitted","Extrapolated"), 
             col = c("gray", "green", "green"),
             lty = c(1,1,2))
    }
    
  }
  
  # ------------------------------------------------------------  
  # NLSB
  
  if (method == "nlsb")
  {
    
    nlsb1 <- fit_nlsb(WSEw_obs)
    
    wbf <- max(cross_section$x)
    wmin <- min(WSEw_obs$w)
    sb.ind <- attributes(nlsb1)$sb.ind
    wsb <- WSEw_obs$w[sb.ind]
    y <- (wbf - wsb)/2
    
    x4 <- seq(wbf/2, (wbf + wmin)/2)
    x5 <- seq((wbf + wmin)/2, wbf-y)
    x6 <- seq(wbf-y, wbf)
    
    b4 <- predict(nlsb1[[1]], newdata = data.frame(w = 2*(x4-wbf/2)))
    b5 <- predict(nlsb1[[1]], newdata = data.frame(w = 2*(x5-wbf/2)))
    # b6 <- predict(nlsb1[[2]], newdata = data.frame(w = 2*(x6-wbf/2)))
    
    ###### This is not working properly
    # z0 <- as.numeric(coef(nlsb1[[1]])[1])# workaround for b6
    # a1 <- as.numeric(coef(nlsb1[[1]])[2]) 
    # s1 <- as.numeric(coef(nlsb1[[1]])[3])
    # a2 <- as.numeric(coef(nlsb1[[2]])[1]) 
    # s2 <- as.numeric(coef(nlsb1[[2]])[2])
    # 
    # predict_nlsb2 <- function(z0, a1, s1, a2, s2, x, offset = wsb)
    # {
    #   pred <- z0 + a1*wsb^s1 + a2*(x-offset)^s2
    # }
    # b6 <- predict_nlsb2(z0, a1, s1, a2, s2, x = 2*(x6-wbf/2))
    #######
    
    # Instead, do this:
    # b5 <- predict(nlsb1[[1]])
    b6 <- predict(nlsb1[[2]])
    # b5 <- approx(b5, n = length(x5))$y
    b6 <- approx(b6, n = length(x6))$y
    
    x1 <- seq(0, y)
    x2 <- seq(y, (wbf - wmin)/2)
    x3 <- seq((wbf - wmin)/2, wbf/2)
    
    b1 <- rev(b6)
    b2 <- rev(b5)
    b3 <- rev(b4)
    
    lines(x1, b1, col = "blue", lty = 1, lwd = 3)
    lines(x2, b2, col = "blue", lty = 1, lwd = 3)
    lines(x3, b3, col = "blue", lty = 2, lwd = 3)
    lines(x4, b4, col = "blue", lty = 2, lwd = 3)
    lines(x5, b5, col = "blue", lty = 1, lwd = 3)
    lines(x6, b6, col = "blue", lty = 1, lwd = 3)
    
    if (showLegend)
    {
      legend("top", legend = c("True","Fitted","Extrapolated"), 
             col = c("gray", "blue", "blue"),
             lty = c(1,1,2))
    }
    
  }
  
}

# Make figure for poster
# par(mfrow = c(2,2), oma = c(2,5,2,2), mar = c(4,5,4,4))
# show_geometry(WSEw = rWSEw[[r]], 
#               WSEw_obs = WSEw_obs1, 
#               cross_section = cross_section,
#               method = "l", ylim = c(134, 145), 
#               showLegend = FALSE, lwd = 3,
#               main = "Linear", 
#               cex.axis = 2, cex.main = 2, cex.lab =2 )
# legend("topleft", legend = c("True","Fitted"), 
#        lty = c(1,1), lwd = c(3,3), 
#        col = c("gray","red"),
#        cex = 1.8, horiz = TRUE)
# show_geometry(WSEw = rWSEw[[r]], 
#               WSEw_obs = WSEw_obs1, 
#               cross_section = cross_section,
#               method = "sb", ylim = c(134, 145), 
#               showLegend = FALSE, lwd = 3,
#               main = "Slope break", 
#               cex.axis = 2, cex.main = 2, cex.lab =2 )
# show_geometry(WSEw = rWSEw[[r]], 
#               WSEw_obs = WSEw_obs1, 
#               cross_section = cross_section,
#               method = "nl", ylim = c(134, 145), 
#               showLegend = FALSE, lwd = 3,
#               main = "Nonlinear", 
#               cex.axis = 2, cex.main = 2, cex.lab =2 )
# show_geometry(WSEw = rWSEw[[r]], 
#               WSEw_obs = WSEw_obs1, 
#               cross_section = cross_section,
#               method = "nlsb", ylim = c(134, 145), 
#               showLegend = FALSE, lwd = 3,
#               main = "Nonlinear slope break", 
#               cex.axis = 2, cex.main = 2, cex.lab =2)

# ------------------------------------------------------------------------------------------------------------------------
# Used a series of if statements, but could also use switch...
# switch(method,
#        "l" = print(1),
#        "sb" = print(2),
#        "nl" = print(3),
#        "nlsb" = print(4)
#        )