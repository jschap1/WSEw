# Find initial guesses
#'
#'Finds initial guesses for the nonlinear model
#' @export
#' @param WSEw_obs data frame containing observations
#' @param extrapolate boolean, if FALSE, uses min observated WSE as initial guess for z0
#' @example init.guess <- make_init_guess(WSEw_obs1)

make_init_guess <- function(WSEw_obs, extrapolate = FALSE)
{
  
  f1 <- lm(WSE ~ I(w^1), data = WSEw_obs) # s = 1
  f2 <- lm(WSE ~ I(w^2), data = WSEw_obs) # s = 2
  f3 <- lm(WSE ~ I(w^3), data = WSEw_obs) # s = 3
  f4 <- lm(WSE ~ I(w^4), data = WSEw_obs) # s = 4
  f5 <- lm(WSE ~ I(w^5), data = WSEw_obs) # s = 5
  f6 <- lm(WSE ~ I(w^6), data = WSEw_obs) # s = 6
  
  s.init <- c(1,2,3,4,5,6)
  
  if (extrapolate)
  {
    z0.init <- as.numeric(c(coef(f1)[1], coef(f2)[1], coef(f3)[1], coef(f4)[1], coef(f5)[1], coef(f6)[1]))
  } else
  {
    z0.init <- rep(min(WSEw_obs$WSE), length(s.init))
  }
  
  a.init <- as.numeric(c(coef(f1)[2], coef(f2)[2], coef(f3)[2], coef(f4)[2], coef(f5)[2], coef(f6)[2]))
  
  init.guess <- data.frame(z0 = z0.init, 
                           a = a.init, 
                           s = s.init) # initial guesses
  
  return(init.guess)
  
}
