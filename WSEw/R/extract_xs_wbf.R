#' Extract cross sections and bankfull width
#'
#' Extracts cross sections and bankfull width
#' @param cross_section output from bisect_line_segments()
#' @param depth bathymetry raster
#' @example result <- extract_xs_wbf(cross_section[39:46], depth)
#' @export
#' @details This is the same as extract_wbf, except that it also outputs the transects
#' and has some error catches.
#' Extracts data underlying the lines in cross_section.
#' Uses first and last overlapped pixels to estimate bankfull width from 
#' perpendicular bisectors to the river centerline.

extract_xs_wbf <- function(cross_section, depth)
{
  n.xs <- length(cross_section)
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
  
  result <- list(wbf = wbf, main_channel = main_channel)
  
  return(result)
}
