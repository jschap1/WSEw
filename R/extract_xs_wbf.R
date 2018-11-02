#' Extract cross sections and bankfull width (parallel)
#'
#' Extracts cross sections and bankfull width in parallel
#' @param cross_section output from bisect_line_segments()
#' @param depth bathymetry raster
#' @param hpc flag for "high performance computing"
#' @param h minimum width (in pixels) to count a channel
#' @examples wbfs_and_mc <- extract_xs_wbf(cross_section, depth, hpc = TRUE)
#' wbfs_and_mc <- extract_xs_wbf(cross_section[1:10], depth, hpc = TRUE)
#' @export
#' @details This is the same as extract_wbf, except that it also outputs the transects
#' and has some error catches.
#' Extracts data underlying the lines in cross_section.
#' Uses first and last overlapped pixels to estimate bankfull width from 
#' perpendicular bisectors to the river centerline.
#' For Pool 5, I got an error about raster I/O. I think this is due to a corrupt raster file.
#' @importFrom foreach foreach
#' @importFrom doMC registerDoMC

extract_xs_wbf <- function(cross_section, depth, hpc = TRUE, h = 20)
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
    
    res1 <- get_main_channel(transects[[1]][,2], return_all = TRUE)
    pixel_widths <- res1$pixel_widths
    
    n.channels <- count_separate_channels_2(transects[[1]][,2], min_width = h)
    if (n.channels > 1)
    {
      wbf_and_mc <- list(NA, NA, NA, n.channels, pixel_widths)
      return(wbf_and_mc)
    }
    
    mc <- get_main_channel(transects[[1]][,2], return_all = TRUE)
    main_channel <- mc$main_channel
    
    first_cell <- transects[[1]][mc$first_ind,1]
    last_cell <- transects[[1]][mc$last_ind,1]
    
    p1 <- xyFromCell(depth, first_cell)
    p2 <- xyFromCell(depth, last_cell)
    
    wbf <- dist(rbind(p1,p2)) 
    
    xy <- rbind(p1,p2)
    xy.sp <- SpatialPoints(xy)
    xs_locations <- SpatialLines(list(Lines(Line(xy.sp), ID = x)))
    
    wbf_and_mc <- list(wbf, main_channel, xs_locations, n.channels, pixel_widths)
    print(paste("Extracted data and estimated wbf for cross section", x, "of", n.xs))
    return(wbf_and_mc)
  }
  
  if (hpc)
  {
    ncores <- detectCores()
    registerDoMC(cores = ncores - 1)
    # registerDoMC(cores = 6)
    print("Set up for parallel processing")
    
    wbfs_and_mc <- foreach(x = 1:n.xs) %dopar% {extract_xs_wbf_par(x)}
    
    # post-process/reformat
    wbf <- vector(length = n.xs) 
    num.channels <- vector(length = n.xs) 
    main_channel <- vector(length = n.xs, "list")
    xs_locations <- vector(length = n.xs, "list")
    pixel_widths <- vector(length = n.xs, "list")
    for (x in 1:n.xs)
    {
      wbf[x] <- wbfs_and_mc[[x]][[1]]
      main_channel[[x]] <- wbfs_and_mc[[x]][[2]]
      xs_locations[[x]] <- wbfs_and_mc[[x]][[3]]
      num.channels[x] <- wbfs_and_mc[[x]][[4]]
      pixel_widths[[x]] <- wbfs_and_mc[[x]][[5]]
    }
    wbfs_and_mc <- list(wbf = wbf, 
                        main_channel = main_channel, 
                        xs_locs = xs_locations, 
                        n.channels = num.channels, 
                        pixel_widths = pixel_widths)
    
  } else # if not hpc
  {
    
    wbf <- vector(length = n.xs)
    n.channels <- vector(length = n.xs)
    transects <- vector(length = n.xs, "list")
    main_channel <- vector(length = n.xs, "list")
    xs_locations <- vector(length = n.xs, "list")
    pixel_widths <- vector(length = n.xs, "list")
    
    for (x in 1:n.xs)
    {
      
      # Output a transect as a KML file
      # xs.sldf <- SpatialLinesDataFrame(cross_section[x], data = data.frame(ID = 1), match.ID = FALSE)
      # xs.sldf.wgs <- spTransform(xs.sldf, crs("+init=epsg:4326"))
      # writeOGR(xs.sldf.wgs, "./Outputs/Cross_Sections/xs_test.kml", layer="reach", driver="KML", overwrite= TRUE)
      
      # try bc extract will fail if  cross_section does not overlap with the extent of depth
      try(transects[[x]] <- extract(depth, cross_section[x], along = TRUE, cellnumbers = TRUE))
      
      # Plot a transect
      # dd <- refWSE[i]*0.3048 - transects[[x]][[1]][,2]
      # plot(dd, type = "l", main = "bed vs x (m)")
      
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
      
      res1 <- get_main_channel(transects[[x]][[1]][,2], return_all = TRUE)
      pixel_widths[[x]] <- res1$pixel_widths
      
      n.channels[x] <- count_separate_channels_2(transects[[x]][[1]][,2], min_width = h)
      if (n.channels > 1)
      {
        wbf[x] <- NA
        main_channel[[x]] <- NA
        xs_locations[[x]] <- NA
        next
      }
      
      mc <- get_main_channel(transects[[x]][[1]][,2], return_all = TRUE)
      main_channel[[x]] <- mc$main_channel
      
      first_cell <- transects[[x]][[1]][mc$first_ind,1]
      last_cell <- transects[[x]][[1]][mc$last_ind,1]
      
      p1 <- xyFromCell(depth, first_cell)
      p2 <- xyFromCell(depth, last_cell)
      
      xy <- rbind(p1,p2)
      xy.sp <- SpatialPoints(xy)
      xs_locations[[x]] <- SpatialLines(list(Lines(Line(xy.sp), ID = x)))
      
      wbf[x] <- dist(rbind(p1,p2)) 
      # wbf[x] <- LineLength(spl@lines[[1]]@Lines[[1]]) # alternative, gives same result
      # estimate of bankfull width, ought to be subject to at most two 
      # pixels worth of error (though it can actually have more error than that)
      
      print(paste("Extracted data and estimated wbf for cross section", x, "of", n.xs))
    }
    
    wbfs_and_mc <- list(wbf = wbf, 
                        main_channel = main_channel, 
                        xs_locs = xs_locations, 
                        n.channels = n.channels, 
                        pixel_widths = pixel_widths)
    
  }
  
  print("Elapsed time is:")
  print(Sys.time() - begin.time)
  return(wbfs_and_mc)
}

# # Find out what are typical widths of side channels
# pixel_widths <- vector(length = n.xs, "list")
# for (x in 1:n.xs)
# {
#   transects[[x]] <- extract(depth, cross_section[x], along = TRUE, cellnumbers = TRUE, progress = "text")
#   if (is.null(transects[[x]][[1]][,2]))
#   {
#     next    
#   }
#   res1 <- get_main_channel(transects[[x]][[1]][,2], return_all = TRUE)
#   pixel_widths[[x]] <- res1$pixel_widths
# }



