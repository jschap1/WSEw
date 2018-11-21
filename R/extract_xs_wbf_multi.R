#' Extract cross sections and bankfull width (multiple channel, parallel)
#'
#' Extracts cross sections and bankfull width in parallel. Sums widths across multichannel transects.
#' @param cross_section output from bisect_line_segments()
#' @param depth bathymetry raster
#' @param hpc flag for "high performance computing." Always TRUE.
#' @param h minimum width (in pixels) to count a channel
#' @export
#' @details See extract_xs_wbf.
#' @importFrom foreach foreach
#' @importFrom doMC registerDoMC

# Cross section x = 2450 has two channels in Pool 21. Use it as a test case.

extract_xs_wbf_multi <- function(cross_section, depth, hpc = TRUE, h = 20)
{
  
  begin.time <- Sys.time()
  print(paste("Start time is ", begin.time))
  
  n.xs <- length(cross_section)
  
  extract_xs_wbf_par <- function(x) # function for evaluating in parallel.
  {
    
    try(transects <- extract(depth, cross_section[x], along = TRUE, cellnumbers = TRUE))
    
    if (!exists("transects")) # in case the cross section is outside the extent of depth (quick fix)
    {
      transects <- NULL
    }
    
    if (is.null(transects[[1]])) # skip empty transects
    {
      print(paste("skipping NULL", x))
      wbf_and_mc <- list(NA, NA, NA, NA, NA)
      return(wbf_and_mc)
    }
    
    if (all(is.na(transects[[1]][,2])))
    {
      print(paste("skipping NA", x))
      wbf_and_mc <- list(NA, NA, NA, NA, NA)
      return(wbf_and_mc)
    }
    
    mc <- get_main_channel(transects[[1]][,2], return_all = TRUE)
    main_channel <- mc$main_channel
    pixel_widths <- mc$pixel_widths
    ch_start_ind <- mc$ind1 # first index of each channel
    ch_end_ind <- mc$ind2 # last index of each channel
    n.channels <- length(ch_start_ind)
    all_channels <- mc$all_channels
  
    # Calculate bankfull width
    # 
    # Calculates the bankfull width of the kth channel of the transect
    # Sums bankfull width across multiple channels (wider than 100 m)
    calc_wbfk <- function(transect, k)
    {
      first_cell <- transects[[1]][mc$ind1[k],1]
      last_cell <- transects[[1]][mc$ind2[k],1]
      p1 <- xyFromCell(depth, first_cell)
      p2 <- xyFromCell(depth, last_cell)
      wbf <- dist(rbind(p1,p2))
      if (wbf <= 100) # bc SWOT can't see narrow rivers
      {
        wbf <- 0
      }
      return(wbf)
    }

    wbf <- 0
    for (k in 1:n.channels)
    {
      wbf <- wbf + calc_wbfk(transects, k)
    }
    
    # make a polyline across the whole cross section (all channels)
    first_cell <- transects[[1]][mc$ind1[1],1] 
    last_cell <- transects[[1]][mc$ind2[n.channels],1]
    p1 <- xyFromCell(depth, first_cell)
    p2 <- xyFromCell(depth, last_cell)
    xy <- rbind(p1,p2)
    xy.sp <- SpatialPoints(xy)
    xs_locations <- SpatialLines(list(Lines(Line(xy.sp), ID = x)))

    wbf_and_mc <- list(wbf, all_channels, xs_locations, n.channels, pixel_widths)
    print(paste("Extracted data and estimated wbf for cross section", x, "of", n.xs))
    return(wbf_and_mc)
  }
  
  ncores <- detectCores()
  registerDoMC(cores = ncores - 1)
  # ncores <- 6
  registerDoMC(cores = ncores)
  print(paste("Set up for parallel processing with", ncores, "cores"))
  
  wbfs_and_mc <- foreach(x = 1:n.xs) %dopar% {extract_xs_wbf_par(x)}
  
  # post-process/reformat
  wbf <- vector(length = n.xs) 
  num.channels <- vector(length = n.xs) 
  all_channels <- vector(length = n.xs, "list")
  xs_locations <- vector(length = n.xs, "list")
  pixel_widths <- vector(length = n.xs, "list")
  for (x in 1:n.xs)
  {
    wbf[x] <- wbfs_and_mc[[x]][[1]]
    all_channels[[x]] <- wbfs_and_mc[[x]][[2]]
    xs_locations[[x]] <- wbfs_and_mc[[x]][[3]]
    num.channels[x] <- wbfs_and_mc[[x]][[4]]
    pixel_widths[[x]] <- wbfs_and_mc[[x]][[5]]
  }
  wbfs_and_mc <- list(wbf = wbf, 
                      all_channels = all_channels, 
                      xs_locs = xs_locations, 
                      n.channels = num.channels, 
                      pixel_widths = pixel_widths)

  print("Elapsed time is:")
  print(Sys.time() - begin.time)
  return(wbfs_and_mc)
}
