#' Calculate mean cross section (main)
#' 
#' @param cross_sections output from auto_transects
#' @param reach_length desired averaging length (m)
#' @param section_length distance between cross sections (m)
#' @return cross_section_avg reach average cross sections
#' @export
#' @example cross_sections_avg <- calc_mean_cross_section(cross_sections, reach_length = 10e3, section_length = 1e3)

calc_mean_cross_section <- function(cross_sections, reach_length, section_length)
{
  
  xs.res <- resample_xs(cross_sections, n = 300) 
  # 300 is arbitrary, just want something fairly large for a small discretization
  
  n.xs <- length(cross_sections$x)
  
  # Get beginning and ending indices of each reach
  ind <- get_start_end_ind(n.xs, reach_length, section_length)
  start.ind <- ind$start.ind
  end.ind <- ind$end.ind
  nr <- ind$nr
  
  na.ind <- get_na_ind(cross_sections)
  
  xs.avg <- vector(length = nr, "list")
  fract <- vector(length = nr)
  for (r in 1:nr)
  {
    xs.avg[[r]] <- calc_mean_cross_section_sub(xs.res[start.ind[r]:end.ind[r],,])
    missing.xs <- sum(na.ind[start.ind[r]:end.ind[r]])
    total.xs <- end.ind[r] - start.ind[r] + 1
    fract[r] <- (total.xs - missing.xs)/total.xs
  }
  
  x <- vector("list", length = nr)
  b <- vector("list", length = nr)
  d <- vector("list", length = nr)
  for (r in 1:nr)
  {
    x[[r]] <- xs.avg[[r]]$x
    b[[r]] <- xs.avg[[r]]$b
    d[[r]] <- xs.avg[[r]]$d
  }

  cross_sections_avg <- list(x = x, b = b, d = d, fract = fract)
  return(cross_sections_avg)
  
}

# ------------------------------------------------------------------------------------------------
#' Calculate Mean Cross Section (subroutine)
#' 
#' A subroutine of calc_mean_cross_section
#' @param xs.res cross section geometry that you wish to average, 
#' resampled to the same length 
#' @export

calc_mean_cross_section_sub <- function(xs.res)
{
  mean.x <- apply(xs.res[,,1], 2, mean, na.rm = TRUE)
  mean.b <- apply(xs.res[,,2], 2, mean, na.rm = TRUE)
  mean.d <- apply(xs.res[,,3], 2, mean, na.rm = TRUE)
  xs.avg <- data.frame(x = mean.x, b = mean.b, d = mean.d)
  return(xs.avg)
}

# ------------------------------------------------------------------------------------------------
#' Resample Cross Section
#' 
#' Uses linear interpolation to resample cross sections with different resolutions
#'to the same resolution.
#' @param cross_sections list of cross section data  
#' @param n resolution of the resample cross sections
#' @examples
#' xs.res <- resample_xs(cross_sections, n)
#' # check that it worked
#' x <- 100
#' plot(cross_sections$x[[x]], cross_sections$b[[x]])
#' lines(xs.res[x,,1], xs.res[x,,2], col = "red")
#' @export

resample_xs <- function(cross_sections, n)
{
  
  na.ind <- get_na_ind(cross_sections)
  
  n.xs <- length(cross_sections$x)
  
  # new data structure where all cross sections have the same number of points
  # data frame would be a better output format than array bc names
  xs.res <- array(dim = c(n.xs, n, 3)) # cross sections by numdata by (x,b,d)
  
  for(x in 1:n.xs) # for each cross section
  {
    
    if ((na.ind[x]) | cross_sections$channel.pix[[x]] <= 30) # really, should restrict in terms of wbf, but this is approximate 
    {
      
      xs.res[x,,] <- cbind(x = rep(NA, n),  
                           b = rep(NA, n), 
                           d = rep(NA, n))
      
    } else
    {
      
      xs <- data.frame(x = cross_sections$x[[x]], b = cross_sections$b[[x]]) # get coordinates of cross section
      d.coord <- cross_sections$d[[x]]
      xs.res[x,,] <- cbind(x = approx(xs, n = n)$x,  # linearly interpolate to the new length n
                           b = approx(xs, n = n)$y, 
                           d = approx(d.coord, n = n)$y)
      
    }

  }
  return(xs.res)
  
}

# ------------------------------------------------------------------------------------------------
#' Get start and end indices of reaches
#' 
#' @param n.xs output from auto_transects
#' @param reach_length desired length of reach (m)
#' @param section_length distance between cross sections (m)
#' @export
#' @example get_start_end_ind(n.xs = 24729, 10e3, 5)

get_start_end_ind <- function(n.xs, reach_length, section_length)
{
  
  # number of full reaches, plus potentially a partial reach that is not the full reach_length long
  n.xs.in.r <- reach_length/section_length # number of cross sections in a reach
  
  if (n.xs %% n.xs.in.r == 0) # if it divides evenly
  {
    print("There are the perfect number of cross sections for division! Lucky you.")
    nr <- n.xs/n.xs.in.r
    start.ind <- seq(1, nr*n.xs.in.r, by = n.xs.in.r)
    end.ind <- start.ind + n.xs.in.r - 1
  } else # if it does not divide evenly
  {
    print("The last cross section is averaged over a smaller distance than you specified")
    nr <- floor(n.xs/n.xs.in.r) + 1 
    print(paste("The averaging distance for cross_section_avg", nr, "is", section_length*(n.xs %% n.xs.in.r), "m"))
    start.ind <- seq(1, n.xs.in.r*(nr-1), by = n.xs.in.r)
    end.ind <- start.ind + n.xs.in.r - 1
    start.ind[nr] <- end.ind[nr-1] + 1
    end.ind[nr] <- n.xs
  }
  
  inds <- list(start.ind = start.ind, end.ind = end.ind, nr = nr)
  return(inds)
  
}