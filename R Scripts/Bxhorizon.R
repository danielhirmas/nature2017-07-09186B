# ---------------------------------------------------------------------------
# Regress the climate data from the stations against an average
# residual effective porosity (x = spatial model) assigned to each station 
# for just the B horizons (aquic moisture regime samples removed)
# ---------------------------------------------------------------------------
#  Finds and calculates the mean EP and residual EP for the closest station
rerun <- T # set to false if ep points to the nearest station are already 
#   found
if (rerun==T) {
  # Checks to see if the climate data.frame already contains these
  # variables and if so, removes them
  climatehold <- which(names(climate)=="Bresepx")
  if (length(climatehold)!=0) {
    climate <- climate[,-climatehold] # removes the offending columns
  }
  #
  # Load the package 'fields' (See note in the preamble above)
  require(fields)
  #
  # Create a distance matrix to find the minimum distance for the 
  #   B horizons
  cl <- as.matrix(data.frame(lon=climate$long,lat=climate$lat))
  da <- as.matrix(data.frame(lon=dBx$long,lat=dBx$lat))
  distmat <- rdist.earth(cl,da) # rows are the row index for cl; cols are the
  #   row index for da
  #
  # Run a loop to record the nearest station to each climate station
  stmat <- matrix(NA,nrow=nrow(distmat),ncol=ncol(distmat))
  for (j in 1:ncol(distmat)) {
    stmat[which.min(distmat[,j]),j] <- 1
  }
  #
  # Run a loop to calculate mean residual effective porosity for the 
  #   B horizons
  steff <- NULL
  for (i in 1:nrow(stmat)) {
    hold <- which(stmat[i,]==1)
    if (length(hold)>0) {
      steff <- c(steff,mean(dBx$RESEPx[hold]))
    } else {
      steff <- c(steff,NA)
    }
  }
  climate <- data.frame(climate,Bresepx=steff)
}
# ---------------------------------------------------------------------------