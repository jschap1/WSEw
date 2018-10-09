#' Get Main Channel
#' 
#' Uses the transect to get the main channel from a transect potentially covering multiple channels. 
#' Also provides other useful info, like width to banks.
#' @param transect depth values along a transect
#' @keywords transects, hydraulics, cross sections
#' @export
#' @examples 
#' get_main_channel(transect)

get_main_channel <- function(transect)
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
  width <- max(summary_table) # width of the main channel
  
  widest.ch <- as.numeric(which.max(summary_table)) # channel number that is the widest
  
  wide.ch.ind <- which(channel == widest.ch) # indices of the widest channel
  
  main_channel <- transect[wide.ch.ind]
  
  return(main_channel)
}
