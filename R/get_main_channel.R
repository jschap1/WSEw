#' Get Main Channel
#' 
#' Uses the transect to get the main channel from a transect potentially covering multiple channels. 
#' Also provides other useful info, like width to banks.
#' @param transect depth values along a transect
#' @param return_all use it to return the indices of transects where the main channel begins and ends
#' @keywords transects, hydraulics, cross sections
#' @export
#' @examples 
#' get_main_channel(transect)
#' @details This function is not perfect!
#' Sometimes the widest channel may not be the main channel.
#' Can we discriminate based on how deep the channel is?
#' With the current method, there is going to be some error 
#' associated with misclassification of the main channel.
#' plot(transect[wide.ch.ind], main = "2", type= "l")
#' plot(transect[wide.ch.ind1],  main = "1", type= "l")
#' summary(transect[wide.ch.ind])
#' summary(transect[wide.ch.ind1])

get_main_channel <- function(transect, return_all = FALSE)
{

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
  
  summary_table <- table(channel) # summarizes the entries in channel
  width <- max(summary_table) # width of the main channel (pixels)
  
  widest.ch <- as.numeric(which.max(summary_table)) # channel number that is the widest
  wide.ch.ind <- which(channel == widest.ch) # indices of the widest channel
  
  # Get start and end indices of all channels
  n.channels <- length(summary_table)
  start_ind <- vector(length = n.channels)
  end_ind <- vector(length = n.channels)
  for (k in 1:n.channels)
  {
    ch_ind <- which(channel == k)
    start_ind[k] <- ch_ind[1]
    end_ind[k] <- ch_ind[length(ch_ind)]
  }
  
  main_channel <- transect[wide.ch.ind]
  
  if (return_all)
  {
    # added in start_ind and end_ind later (11/19/2018)
    # they indicate the beginnings and ends of each channel, not just the widest one, 
    # which is what first_ind and last_ind do.
    main_channel <- list(main_channel = main_channel, 
                         first_ind = wide.ch.ind[1], 
                         last_ind = wide.ch.ind[length(wide.ch.ind)], 
                         pixel_widths = summary_table, 
                         ind1 = start_ind, 
                         ind2 = end_ind, 
                         all_channels = transect)
  }
  
  return(main_channel)
}


