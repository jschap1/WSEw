# Helper functions for auto_transects workflow
# June 22, 2018 JRS

# ------------------------------------------------------------------------------------

get_main_channel <- function(transect)
{
  # Uses the transect to get the main channel
  # Also provides other useful info, like width to banks
  
  # Assign a value to each river channel intersected by the transect
  
  n <- length(transect)
  channel <- vector(length = n)
  
  ch <- 0
  
  # initialize with one go-through
  if (is.na(transect[1]))
  {
    channel[1] <- NA
  } else 
  {
    ch <- ch + 1
    channel[1] <- ch
  }
  
  # loop over the rest of the entries in transect
  for (i in 2:n)
  {
    if (length(transect[i]!=0) & length(transect[i-1]!=0))
    {
      if (!is.na(transect[i]) & is.na(transect[i-1]))
      {
        ch <- ch + 1
      }

      if (is.na(transect[i]))
      {
        channel[i] <- NA
      } else
      {
       channel[i] <- ch
      }
    }
  }
  
  # Find the widest channel
  # channel.f <- factor(channel, levels = c(1,2))
  
  summary_table <- table(channel) # summarizes the entries in channel
  width <- max(summary_table) # width of the main channel
  
  widest.ch <- as.numeric(which.max(summary_table)) # channel number that is the widest
  
  wide.ch.ind <- which(channel == widest.ch) # indices of the widest channel
  
  main_channel <- transect[wide.ch.ind]
  
  return(main_channel)
}

# ------------------------------------------------------------------------------------

find_closest <- function(rpolyline, polyline, i)
{ # finds the closest row of polyline to the row of rpolyline
  # for use in interpolating widths.
  rpolyline[i,] # find the closest row in polyline
  dist1 <- vector(length = dim(polyline)[1])
  for (j in 1:dim(polyline)[1])
  {
    dist1[j] <- sqrt((polyline$x[j] - rpolyline$x[i])^2 + (polyline$y[j] - rpolyline$y[i])^2)
  }
  r <- which.min(dist1)
  return(r)
}

# ------------------------------------------------------------------------------------

bisect_line_segments <- function(rpolyline, projcrs, depth, w)
{
  # Bisects each line segment in the polyline
  
  w <- w/5 # width in pixels, resolution is 5 m, numerator is in meters 
  
  nseg <- dim(rpolyline)[1]-1
  coords.df <- array(dim = c(nseg, 4))
  coords.df <- as.data.frame(coords.df)
  colnames(coords.df) <- c("x1","y1","x2","y2")
  
  lines_obj <- vector("list", length = nseg)
  
  for (seg in 1:nseg)
  {
    theta <- atan((rpolyline$y[seg+1] - rpolyline$y[seg])/(rpolyline$x[seg+1] - rpolyline$x[seg]))
    if(is.na(theta)) {theta <- pi/2} # in case of vertical segment
    theta_perp <- theta + pi/2
    
    # Draw the bisector
    resolution <- res(depth)[1]
    xb1 <- (rpolyline$x[seg] + rpolyline$x[seg+1])/2
    yb1 <- (rpolyline$y[seg] + rpolyline$y[seg+1])/2
    
    # Extend it out across the river channel
    # (just needs to be at least half as wide as the widest part of the river)
    xb0 <- xb1 - w*resolution*cos(theta_perp)
    xb2 <- xb1 + w*resolution*cos(theta_perp)
    yb0 <- yb1 - w*resolution*sin(theta_perp)
    yb2 <- yb1 + w*resolution*sin(theta_perp)
    
    coords <- cbind(c(xb0, xb2), c(yb0, yb2)) # each line is defined by two points
    coords.df[seg,] <- c(coords[1,1], coords[1,2], coords[2,1], coords[2,2])
    # Outputs a data frame with coordinates for the endpoints of the bisectors of rpolyline
    
    # Make a SpatialLines object for each bisector
    simple_line <- Line(cbind(rbind(coords.df$x1[seg], coords.df$x2[seg]),
                              rbind(coords.df$y1[seg], coords.df$y2[seg])))
    lines_obj[[seg]] <- Lines(list(simple_line), ID = as.character(seg)) # can list these
    
  }
  cross_section <- SpatialLines(lines_obj,proj4string = projcrs)
  return(cross_section)
}

# ------------------------------------------------------------------------------------

get_depth <- function(transect)
{
  lutable <- depth@data@attributes[[1]]
  transect.depth <- vector(length = length(transect))
  for (k in 1:length(transect))
  {
    transect.depth[k] <- as.numeric.factor(lutable$DEPTH_M[match(transect[k],lutable$ID)])
    if (k%%50000 == 0)
    {
      print(k)
    }
  }
  return(transect.depth)
}

# ------------------------------------------------------------------------------------

estimate_widths <- function(rpolyline, resolution = 5, tlength, nseg)
{
  # Uses angle and number of pixels crossed to estimate cross-section width
  # Simplification. Should check if it is OK assumption.
  cross_section_width <- vector(length=nseg)
  for (seg in 1:nseg)
  {
    theta <- atan((rpolyline$y[seg+1] - rpolyline$y[seg])/(rpolyline$x[seg+1] - rpolyline$x[seg]))
    if(is.na(theta)) {theta <- pi/2} # in case of vertical segment
    theta_perp <- theta + pi/2
    cross_section_width[seg] <- tlength[[seg]]*abs(cos(theta_perp))*resolution
  }
  return(cross_section_width)
}

# ------------------------------------------------------------------------------------

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

