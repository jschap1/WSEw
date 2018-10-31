#' Subset cross sections
#' 
#' Subsets the variable cross_sections using an index, ind
#' @param cross_sections from build_cross_sections
#' @param ind index of cross sections to keep
#' @example cross_sections_sub <- subset_cross_sections(cross_sections, ind)

subset_cross_sections <- function(cross_sections, ind)
{
  n <- length(ind)
  sub_xs <- list()
  sub_xs$x <- vector(length = n, "list")
  sub_xs$b <- vector(length = n, "list")
  sub_xs$d <- vector(length = n, "list")
  count <- 1
  for (i in 1:length(ind))
  {
    sub_xs$x[[count]] <- cross_sections$x[[ind[i]]]
    sub_xs$b[[count]] <- cross_sections$b[[ind[i]]]
    sub_xs$d[[count]] <- cross_sections$d[[ind[i]]]
    count <- count + 1
  }
  return(sub_xs)
}

# nbr.ind <- c(3, 4, 5)
# csa.nbr <- subset_cross_sections(cross_sections_avg, nbr.ind)
# nr <- length(csa.nbr)

# ------------------------------------------------------------------------------------------------------------

#' Subset cross sections
#' 
#' Subsets the variable cross_sections using an index, ind
#' @param cross_sections from build_cross_sections
#' @param ind index of cross sections to keep
#' @example rWSEw_sub <- subset_rWSEw(rWSEw, nbr.ind)

subset_rWSEw <- function(rWSEw, ind)
{
  n <- length(ind)
  sub_WSEw <- vector(length = n, "list")
  count <- 1
  for (i in 1:length(ind))
  {
    sub_WSEw[[count]]$w <- rWSEw[[ind[i]]]$w
    sub_WSEw[[count]]$WSE <- rWSEw[[ind[i]]]$WSE
    count <- count + 1
  }
  sub_WSEw <- lapply(sub_WSEw, as.data.frame)
  return(sub_WSEw)
}
