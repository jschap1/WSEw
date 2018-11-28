#' Model one reach
#' 
#' Models one reach, using one or multiple methods
#' @export
#' @param full_expo boolean for whether or not to exclude observations at less than 100 m width
#' @return A0 predictions from each method

model_one_reach <- function(exp_dir, r, expo, 
                            sd_wse, sd_w, 
                            method = "all", 
                            full_expo = FALSE,
                            ...)
{
  
  rWSEw <- readRDS(file.path(exp_dir, "rWSEw.rds"))
  
  if (full_expo)
  {
    WSEw_obs1 <- observe2(WSEw = rWSEw[[r]], exposure = expo)  
  }
  else
  {
    WSEw_obs1 <- observe(WSEw = rWSEw[[r]], exposure = expo, sd_wse = sd_wse, sd_w = sd_w)  
  }
  
  
  plot(WSE~w, WSEw_obs1, xlab = "w (m)", ylab = "WSE (m)", ...)
  
  # True h-w line
  lines(WSE~w, rWSEw[[r]], col = "black", lty = 1, lwd = 0.5)
  
  # plot linear
  if (method == "l" | method == "all")
  {
    lf1 <- fit_linear(WSEw_obs1, h = 5)
    lines(WSEw_obs1$w, predict(lf1), col = "red")
    w.new <- seq(0, min(WSEw_obs1$w), length.out = 50)
    lines(w.new, predict(lf1, newdata = data.frame(w = w.new)), col = "red", lty = 2)
    points(0, predict(lf1, newdata = data.frame(w=0)), col = "red", pch = 19)
    if (method == "l")
    {
      return(lf1)
    }
  }

  # plot slope break
  if (method == "sb" | method == "all")
  {
    sb1 <- fit_slopebreak(WSEw_obs1, multiple_breaks = FALSE, continuity = TRUE)
    nn <- length(WSEw_obs1$w)
    sb1.ind <- attributes(sb1)$sb.ind
    print(expo)
    print(sb1.ind)
    lines(WSEw_obs1$w[1:sb1.ind], predict(sb1[[1]]), col = "orange")
    lines(WSEw_obs1$w[(sb1.ind):nn], predict(sb1[[2]]), col = "orange")
    w.new <- seq(0, min(WSEw_obs1$w[1:sb1.ind]), length.out = 50)
    lines(w.new, predict(sb1[[1]], newdata = data.frame(w = w.new)), col = "orange", lty = 2)
    points(0, predict(sb1[[1]], newdata = data.frame(w=0)), col = "orange", pch = 19)
    if (method == "sb")
    {
      return(sb1)
    }
  }


  # plot nl
  if (method == "nl" | method == "all")
  {
    nl1 <- fit_nonlinear(WSEw_obs1)
    lines(WSEw_obs1$w, predict(nl1), col = "green")
    w.new <- seq(0, min(WSEw_obs1$w), length.out = 50)
    lines(w.new, predict(nl1, newdata = data.frame(w = w.new)), col = "green", lty = 2)
    points(0, predict(nl1, newdata = data.frame(w=0)), col = "green", pch = 19)
    if (method == "nl")
    {
      return(nl1)
    }
  }
  
  # plot nlsb
  if (method == "nlsb" | method == "all")
  {
    nn <- length(WSEw_obs1$w)
    nlsb1 <- fit_nlsb(WSEw_obs1)
    nlsb1.ind <- attributes(nlsb1)$sb.ind
    lines(WSEw_obs1$w[1:nlsb1.ind], predict(nlsb1[[1]]), col = "blue")
    lines(WSEw_obs1$w[(nlsb1.ind):nn], predict(nlsb1[[2]]), col = "blue")
    w.new <- seq(0, min(WSEw_obs1$w[1:nlsb1.ind]), length.out = 50)
    lines(w.new, predict(nlsb1[[1]], newdata = data.frame(w = w.new)), col = "blue", lty = 2)
    points(0, predict(nlsb1[[1]], newdata = data.frame(w=0)), col = "blue", pch = 19)
    if (method == "nlsb")
    {
      return(nlsb1)
    }
  }
  
  # Calculate A0 predictions from each model:
  if (method == "all")
  {
    
    # convert format (NL)
    model1 <- nl1
    nl11 <- list(z0 = as.numeric(coef(model1))[1],
                         a = as.numeric(coef(model1))[2],
                         s = as.numeric(coef(model1))[3],
                         WSEbf = max(fitted(model1)),
                         ef = attributes(model1)$ef)
    
    # convert format (NLSB)
    model1 <- nlsb1
    nlsb11 <- list(z0 = as.numeric(coef(model1[[1]]))[1],
         a1 = as.numeric(coef(model1[[1]]))[2],
         s1 = as.numeric(coef(model1[[1]]))[3],
         a2 = as.numeric(coef(model1[[2]]))[1],
         s2 = as.numeric(coef(model1[[2]]))[2],
         WSEbf = max(fitted(model1[[2]])),
         sb.ind = attributes(model1)$sb.ind,
         ef = attributes(model1)$ef)
    
    A0.l1 <- calc_model_A0(lf1, type = "linear")
    A0.sb1 <- calc_model_A0(sb1, type = "sb")
    A0.nl1 <- calc_model_A0(nl11, type = "nl", WSEw_obs = WSEw_obs1)
    A0.nlsb1 <- calc_model_A0(nlsb11, type = "nlsb", WSEw_obs = WSEw_obs1)
    A0 <- list(l = as.numeric(A0.l1), 
               sb = as.numeric(A0.sb1), 
               nl = as.numeric(A0.nl1), 
               nlsb = as.numeric(A0.nlsb1)
               )
    z0 <- list(l = as.numeric(coef(lf1)[1]), 
               sb = as.numeric(coef(sb1[[1]])[1]), 
               nl = as.numeric(coef(nl1)[1]), 
               nlsb = as.numeric(coef(nlsb1[[1]])[1])
               )
    result <- list(z0 = z0, A0 = A0)
    
  } else
  {
    A0 <- NULL
  }
  
  return(result)
  
}
