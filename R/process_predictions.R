#' Process 
#' 
#' Processes L, SB, NL, and NLSB predictions to a more convenient format
#' @export

process_predictions <- function(pred_lf, pred_sb, pred_nl, pred_nlsb, saveloc)
{
  
  nr <- length(pred_lf)
  n_exp_levels <- dim(pred_lf[[1]]$z0)[1]
  M <- dim(pred_lf[[1]]$z0)[2]
  
  # Initialize
  z0.l <- array(dim = c(nr, n_exp_levels, M))
  z0.sb <- array(dim = c(nr, n_exp_levels, M))
  z0.nl <- array(dim = c(nr, n_exp_levels, M))
  z0.nlsb <- array(dim = c(nr, n_exp_levels, M))
  
  A.l <- array(dim = c(nr, n_exp_levels, M))
  A.sb <- array(dim = c(nr, n_exp_levels, M))
  A.nl <- array(dim = c(nr, n_exp_levels, M))
  A.nlsb <- array(dim = c(nr, n_exp_levels, M))
  
  WP.l <- array(dim = c(nr, n_exp_levels, M))
  WP.sb <- array(dim = c(nr, n_exp_levels, M))
  WP.nl <- array(dim = c(nr, n_exp_levels, M))
  WP.nlsb <- array(dim = c(nr, n_exp_levels, M))
  
  A0.l <- array(dim = c(nr, n_exp_levels, M))
  A0.sb <- array(dim = c(nr, n_exp_levels, M))
  A0.nl <- array(dim = c(nr, n_exp_levels, M))
  A0.nlsb <- array(dim = c(nr, n_exp_levels, M))
  
  ef.l <- array(dim = c(nr, n_exp_levels, M))
  ef.sb <- array(dim = c(nr, n_exp_levels, M))
  ef.nl <- array(dim = c(nr, n_exp_levels, M))
  ef.nlsb <- array(dim = c(nr, n_exp_levels, M))
  
  for (r in 1:nr)
  {
    z0.l[r,,] <- pred_lf[[r]]$z0
    z0.sb[r,,] <- pred_sb[[r]]$z0
    z0.nl[r,,] <- pred_nl[[r]]$z0
    z0.nlsb[r,,] <- pred_nlsb[[r]]$z0
    
    A0.l[r,,] <- pred_lf[[r]]$A0
    A0.sb[r,,] <- pred_sb[[r]]$A0
    A0.nl[r,,] <- pred_nl[[r]]$A0
    A0.nlsb[r,,] <- pred_nlsb[[r]]$A0
    
    A.l[r,,] <- pred_lf[[r]]$A
    A.sb[r,,] <- pred_sb[[r]]$A
    A.nl[r,,] <- pred_nl[[r]]$A
    A.nlsb[r,,] <- pred_nlsb[[r]]$A
    
    WP.l[r,,] <- pred_lf[[r]]$WP
    WP.sb[r,,] <- pred_sb[[r]]$WP
    WP.nl[r,,] <- pred_nl[[r]]$WP
    WP.nlsb[r,,] <- pred_nlsb[[r]]$WP
    
    ef.l[r,,] <- pred_lf[[r]]$ef
    ef.sb[r,,] <- pred_sb[[r]]$ef
    ef.nl[r,,] <- pred_nl[[r]]$ef
    ef.nlsb[r,,] <- pred_nlsb[[r]]$ef
  }
  
  z0name <- file.path(saveloc, "z0_pred.rda")
  Aname <- file.path(saveloc, "A_pred.rda")
  WPname <- file.path(saveloc, "WP_pred.rda")
  A0name <- file.path(saveloc, "A0_pred.rda")
  efname <- file.path(saveloc, "ef_pred.rda")
  
  save(z0.l, z0.sb, z0.nl, z0.nlsb, file = z0name)
  save(A.l, A.sb, A.nl, A.nlsb, file = Aname)
  save(WP.l, WP.sb, WP.nl, WP.nlsb, file = WPname)
  save(A0.l, A0.sb, A0.nl, A0.nlsb, file = A0name)
  save(ef.l, ef.sb, ef.nl, ef.nlsb, file = efname)
  
  print(paste("Saved", z0name))
  print(paste("Saved", Aname))
  print(paste("Saved", WPname))
  print(paste("Saved", A0name))
  print(paste("Saved", efname))
  
  z0 <- list(l = z0.l, sb = z0.sb, nl = z0.nl, nlsb = z0.nlsb)
  A <- list(l = A.l, sb = A.sb, nl = A.nl, nlsb = A.nlsb)
  WP <- list(l = WP.l, sb = WP.sb, nl = WP.nl, nlsb = WP.nlsb)
  A0 <- list(l = A0.l, sb = A0.sb, nl = A0.nl, nlsb = A0.nlsb)
  ef <- list(l = ef.l, sb = ef.sb, nl = ef.nl, nlsb = ef.nlsb)
  
  # Compute slope via finite difference
  s0.l <- apply(z0.l, c(2,3), diff)
  s0.sb <- apply(z0.sb, c(2,3), diff)
  s0.nl <- apply(z0.nl, c(2,3), diff)
  s0.nlsb <- apply(z0.nlsb, c(2,3), diff)
  s0 <- list(l = s0.l, sb = s0.sb, nl = s0.nl, nlsb = s0.nlsb)
  s0name <- file.path(saveloc, "s0_pred.rda")
  save(s0.l, s0.sb, s0.nl, s0.nlsb, file = s0name)
  print(paste("Saved", s0name))
  
  return(list(z0 = z0, A = A, WP = WP, A0 = A0, ef = ef, s0 = s0, 
              fnames = c(z0name, Aname, WPname, A0name, efname, s0name)))
  
}
  