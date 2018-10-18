#' Load Pepsi data from text file
#' 
#' @details Saves two versions of the h-w data: one is sorted from low to high width, 
#' and the other is left in the original simulated order of day 1 to nt.
#' @export
#' @importFrom tools file_path_sans_ext
#' Loads width, height, and flow area data from Pepsi Challenge outputs into R
#' Data should be saved as text file, e.g. generated using preprocess_pepsi_data.m

load_pepsi_data_from_text <- function(data.dir, save.dir)
{
  WSEnames <- list.files(data.dir, glob2rx("*_WSE.txt"))
  wnames <- list.files(data.dir, glob2rx("*_w.txt"))
  Anames <- list.files(data.dir, glob2rx("*_A.txt"))
  
  for (i in 1:length(wnames))
  {
    
    h <- t(read.table(file.path(data.dir, WSEnames[i]), sep = ',', header = FALSE, row.names = NULL))
    w <- t(read.table(file.path(data.dir, wnames[i]), sep = ',', header = FALSE, row.names = NULL))
    A <- t(read.table(file.path(data.dir, Anames[i]), sep = ',', header = FALSE, row.names = NULL))
    
    nr <- dim(h)[1] # number of reaches
    
    # Reformat into rWSEw form
    rWSEw <- vector(length = nr, "list")
    for (r in 1:nr)
    {
      wr <- w[r,] # Sort them from low to high WSE
      hr <- h[r,]
      rWSEw.unsorted[[r]] <- data.frame(WSE = hr, w = wr)
      rWSEw.sorted <- sortWSEw(rWSEw.unsorted)
    }
    
    nametmp <- file_path_sans_ext(WSEnames[i]) # get river name
    rivname <- strsplit(nametmp, split = '_WSE')[[1]]
    saveRDS(rWSEw.sorted, file = file.path(save.dir, paste0(rivname, '_rWSEw_sorted.rds')))
    saveRDS(rWSEw.unsorted, file = file.path(save.dir, paste0(rivname, '_rWSEw_unsorted.rds')))
    saveRDS(A, file = file.path(save.dir, paste0(rivname, '_A.rds')))
    print(paste("Saved files to", save.dir))
    
  }
}

# ------------------------------------------------------------------------------------------------
#' Sort WSEw
#' 
#' Sorts the WSEw data in increasing order of flow width
#' @export
#' @example WSEw.sorted <- sortWSEw(WSEw.unsorted)
sortWSEw <- function(WSEw)
{
  
  nr <- length(WSEw)
  rWSEw <- vector(length = nr, "list") 
  for (r in 1:nr)
  {
    hr <- WSEw[[r]]$WSE
    wr <- WSEw[[r]]$w
    df <- data.frame(WSE = hr[order(wr)], w = wr[order(wr)])
    rWSEw[[r]] <- df
  }
  
  return(rWSEw)
  
}
