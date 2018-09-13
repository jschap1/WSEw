#' Test for Gaussianity
#' 
#' Uses the Shapiro-Wilks test to test for Gaussianity
#' @param x data
#' @param alpha significance level
#' @return z 1 if data are Gaussian, 0 otherwise
#' @export

is_normal <- function(x)
{
  if (sum(!is.na(x))<3) # check that there are data
  {
    z <- NA
    return(z)
  } else
  {
    t <- shapiro.test(x)
    if (t$p.value > 0.05) # Gaussian
    {
      z <- 1
    } else # not Gaussian
    {
      z <- 0
    }
    return(z)
  }
}

# And here are some visual tests ----------------------------------------------------------------------------
# summary(z0.l[2,19,])
# hist(z0.l[5,19,], "fd")
# 
# # visual tests
# library(ggplot2)
# library(ggpubr)
# ggdensity(z0.l[5,19,], main = "Density plot")
# ggqqplot(z0.l[5,19,])
# 
# # Shapiro-Wilk's method
# shapiro.test(z0.l[5,19,]) # if p>0.05, then we cannot reject the null hypothesis that the data are normal
