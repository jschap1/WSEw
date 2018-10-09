# Pepsi Challenge data
#
# Oct. 9, 2018
# Expanding the study area
# Read the h-w data using Mike's Matlab code, then save as a text file to load into R.

library(tools)

setwd("/Users/jschap/Documents/Research/SWOTBATH")

data.dir <- "./Data/Pepsi_Challenge_I/WSEw_Pepsi1"

WSEnames <- list.files(data.dir, glob2rx("*_WSE.txt"))
wnames <- list.files(data.dir, glob2rx("*_w.txt"))

for (i in 1:length(wnames))
{
  
  h <- t(read.table(file.path(data.dir, WSEnames[i]), sep = ',', header = FALSE, row.names = NULL))
  w <- t(read.table(file.path(data.dir, wnames[i]), sep = ',', header = FALSE, row.names = NULL))
  
  nr <- dim(h)[2] # number of reaches

  # Reformat into rWSEw form
  rWSEw <- vector(length = nr, "list")
  for (r in 1:nr)
  {
    wr <- w[,r] # Sort them from low to high WSE
    hr <- h[,r]
    df <- data.frame(WSE = hr[order(hr)], w = wr[order(hr)])
    rWSEw[[r]] <- df
  }
  
  nametmp <- file_path_sans_ext(WSEnames[i]) # get river name
  rivname <- strsplit(nametmp, split = '_WSE')[[1]]
  # plot(WSE~w, rWSEw[[1]])
  
  saveRDS(rWSEw, file = paste0(rivname, '_rWSEw.rds'))
  
}

# Now, can fit statistical models to the h-w data saved as the output of this script.
