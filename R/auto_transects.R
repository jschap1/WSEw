#' Auto Transects
#' 
#' Draws transects, finds main channel, and gets depth values
#' @export

auto_transects <- function(section_length, riv, depth, halfwidth, saveflag = FALSE, output_all = FALSE, savedir = NULL)
{

  rivsplit <- splitLines(riv, dist = section_length)
  projcrs <- crs(depth)
  crs(rivsplit) <- projcrs

  # Bisect line segments
  cross_section <- bisect_line_segments2(rivsplit, w = halfwidth, resolution = res(depth)[1], mid = TRUE)

  # Make polyline containing only the ending transects
  n.xs <- length(cross_section)
  ind <- get_start_end_ind(n.xs, reach_length = 10e3, section_length)
  start.ind <- ind$start.ind
  end.ind <- ind$end.ind
  xs.end <- cross_section[end.ind]

  # Plot study area map
  plot(depth, main = "Bathymetry", xlab = "Easting", ylab = "Northing", legend = FALSE)
  lines(riv, col = "blue")
  lines(cross_section, col = "Red", lwd = 0.5) # the segments are numbered from north to south
  lines(xs.end, col = "black", lwd = 2)

  print("Extracting values along transects.")
  print("This can take a LONG time.")
  print("It takes about 12 seconds per cross section per processor")
  res <- extract_xs_wbf(cross_section, depth, hpc = TRUE, h = 20)
  wbf <- res$wbf
  main_channel <- res$main_channel
  xs_locs <- res$xs_locs # a single SpatialLines object with IDs corresponding to XS number
  pixel_widths <- res$pixel_widths
  n.channels <- res$n_channels

  # Combine the SpatialLines objects together
  # Need to remove NA objects prior to combining
  rm.ind <- which(wbf == 0 | is.na(wbf))
  if (any(rm.ind))
  {
    xs_locs <- xs_locs[-rm.ind] # Slot ID keeps track of the cross section
  }
  xs_locs_combined <- do.call(rbind, xs_locs)

  if (saveflag)
  {
    saveRDS(cross_section, file = file.path(savedir, "cross_sections_polyline.rds"))
    saveRDS(xs.end, file = file.path(savedir, "reach_ends.rds"))
    saveRDS(wbf, file = file.path(savedir, "wbf.rds")) # bankfull width
    saveRDS(main_channel, file = file.path(savedir, "main_channel.rds")) # depth of main channel
    saveRDS(xs_locs_combined, file = file.path(savedir, "xs_locations.rds")) # cross section locations
    saveRDS(pixel_widths, file = file.path(savedir, "pixel_widths.rds")) # width of cross sections, in pixels
    saveRDS(n.channels, file = file.path(savedir, "n_channels.rds")) # number of channels crossed by each transect
  }

  if (output_all)
  {
    result <- list(d = main_channel, wbf = wbf)
    return(result)
  } else
  {
    d <- main_channel
    return(d)
  }

}
