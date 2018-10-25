#' Plot height-width model 2
#' 
#' For making movies of the fitted models with increasing channel exposure
#' @param WSEw h-w data
#' @param r reach number
#' @param z0 z0 predictions, must correpond to the type of model
#' @param type type of model. Can be "l", "sb", "nl", or "nlsb"
#' @param k exposure level index
#' @param ... additional arguments to plot()
#' @details Requires that models have already been fitted and saved as .rds files, 
#' and that predicted hydraulic parameters (specifically z0) have already been computed.
#' Must run from the SWOTBATH/home directory.
#' Currently there is no option for no measurement error.
#' Currently only has methods for linear model.
#' @export
#' @examples plot_hw_fits_2(rWSEw, z0.l, r = 1, type = "l", ylim=c(138,143), exp_dir = exp_dir, k = 8)

plot_hw_fits_2 <- function(WSEw, z0, r, type, errorflag = TRUE, exp_dir, k, ...)
{
  
  plot(WSE~w, WSEw[[r]], type = "l", ...)
  
  if (errorflag == TRUE)
  {
    obsname <- paste0("obs/WSEw_obs_r_", r, ".rds")
    WSEw_obs <- readRDS(file.path(exp_dir, obsname))
    print(k)
  } else
  {
    WSEw_obs <- observe(WSEw[[r]], exposure = 1, sd_wse = 0, sd_w = 0)
  }
  
  M <- length(WSEw_obs[[k]])
  
  for (m in 1:M)
  {
    points(WSE~w, WSEw_obs[[k]][[m]], col = "cyan")
  }
  
  if (type == "l")
  {
    modelname <- paste0("lf/lf_", "r_", r, ".rds")
    model <- readRDS(file.path(exp_dir, modelname))
    
    for (m in 1:M)
    {
      lines(WSEw_obs[[k]][[m]]$w, predict(model[[k]][[m]]), col = "red")  
    }
    
    points(rep(0, M), z0.l[r,k,], col = "red", pch = 19, cex = 1)
    
  }
  
  if (type == "sb")
  {
    modelname <- paste0("sb/sb_", "r_", r, ".rds")
    model <- readRDS(file.path(exp_dir, modelname))
  }
  
  
  if (type == "nl")
  {
    modelname <- paste0("nl/nl_", "r_", r, ".rds")
    model <- readRDS(file.path(exp_dir, modelname))
  }
  
  
  if (type == "nlsb")
  {
    modelname <- paste0("nlsb/nlsb_", "r_", r, ".rds")
    model <- readRDS(file.path(exp_dir, modelname))
  }
  
}