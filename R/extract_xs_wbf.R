#' Extract cross sections and bankfull width (parallel)
#'
#' Extracts cross sections and bankfull width in parallel
#' @param cross_section output from bisect_line_segments()
#' @param depth bathymetry raster
#' @param hpc flag for "high performance computing"
#' @examples wbfs_and_mc <- extract_xs_wbf(cross_section, depth, hpc = TRUE)
#' wbfs_and_mc <- extract_xs_wbf(cross_section[1:10], depth, hpc = TRUE)
#' @export
#' @details This is the same as extract_wbf, except that it also outputs the transects
#' and has some error catches.
#' Extracts data underlying the lines in cross_section.
#' Uses first and last overlapped pixels to estimate bankfull width from 
#' perpendicular bisectors to the river centerline.
#' @importFrom foreach foreach
#' @importFrom doMC registerDoMC

extract_xs_wbf <- function(cross_section, depth, hpc = TRUE)
{
  
  begin.time <- Sys.time()
  print(paste("Start time is ", begin.time))
  
  n.xs <- length(cross_section)

  extract_xs_wbf_par <- function(x) # function for evaluating in parallel.
  {
    transects <- extract(depth, cross_section[x], along = TRUE, cellnumbers = TRUE)
    
    if (is.null(transects[[1]])) # skip empty transects
    {
      print(paste("skipping NULL", x))
      wbf_and_mc <- list(NA, NA)
      return(wbf_and_mc)
    }
    
    if (all(is.na(transects[[1]][,2])))
    {
      print(paste("skipping NA", x))
      wbf_and_mc <- list(NA, NA)
      return(wbf_and_mc)
    }
    
    mc <- get_main_channel(transects[[1]][,2], return_index = TRUE)
    main_channel <- mc$main_channel
    
    first_cell <- transects[[1]][mc$first_ind,1]
    last_cell <- transects[[1]][mc$last_ind,1]
    
    p1 <- xyFromCell(depth, first_cell)
    p2 <- xyFromCell(depth, last_cell)
    
    wbf <- dist(rbind(p1,p2)) 
    
    wbf_and_mc <- list(wbf, main_channel)
    print(paste("Extracted data and estimated wbf for cross section", x, "of", n.xs))
    return(wbf_and_mc)
  }
  
  if (hpc)
  {
    ncores <- detectCores()
    registerDoMC(cores = ncores - 1)
    print("Set up for parallel processing")
    
    wbfs_and_mc <- foreach(x = 1:n.xs) %dopar% {extract_xs_wbf_par(x)}
    
    # post-process/reformat
    wbf <- vector(length = n.xs) 
    main_channel <- vector(length = n.xs, "list")
    for (x in 1:n.xs)
    {
      wbf[x] <- wbfs_and_mc[[x]][[1]]
      main_channel[[x]] <- wbfs_and_mc[[x]][[2]]
    }
    wbfs_and_mc <- list(wbf = wbf, main_channel = main_channel)
    
  } else
  {
    
    wbf <- vector(length = n.xs)
    transects <- vector(length = n.xs, "list")
    main_channel <- vector(length = n.xs, "list")
    
    for (x in 1:n.xs)
    {
      transects[[x]] <- extract(depth, cross_section[x], along = TRUE, cellnumbers = TRUE)
      
      if (is.null(transects[[x]][[1]])) # skip empty transects
      {
        print(paste("skipping NULL", x))
        next
      }
      
      if (all(is.na(transects[[x]][[1]][,2])))
      {
        print(paste("skipping NA", x))
        next
      }
      
      mc <- get_main_channel(transects[[x]][[1]][,2], return_index = TRUE)
      main_channel[[x]] <- mc$main_channel
      
      first_cell <- transects[[x]][[1]][mc$first_ind,1]
      last_cell <- transects[[x]][[1]][mc$last_ind,1]
      
      p1 <- xyFromCell(depth, first_cell)
      p2 <- xyFromCell(depth, last_cell)
      
      wbf[x] <- dist(rbind(p1,p2)) 
      # estimate of bankfull width, ought to be subject to at most two 
      # pixels worth of error (though it can actually have more error than that)
      
      print(paste("Extracted data and estimated wbf for cross section", x, "of", n.xs))
    }
    
    wbfs_and_mc <- list(wbf = wbf, main_channel = main_channel)
    
  }
  
  print("Elapsed time is:")
  print(Sys.time() - begin.time)
  return(wbfs_and_mc)
}

