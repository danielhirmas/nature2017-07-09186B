# ---------------------------------------------------------------------------
# Examine the relationship between climate and ep for A horizons (see
# the selection of A horizons below) that have their mid horizon depth within 
# the upper 25 cm
# ---------------------------------------------------------------------------
setwd(fname)
# Select the A horizons
mtlt <- "AHorizonswithMidpointDepthswithin25cm"
mindepth <- 25 # selects horizons that are less than or equal to this depth 
# (cm)
horizonhold <- c("A","A2","A1","A3","A12","A11","Ac1","Ag","A4","Ag1",
                 "Ac2","A21","A13","Ass","Ag2","A1c","A22")
require(fields)
hhold <- NULL
for (i in 1:nrow(dat)) {
  if ((length(grep("A",dat$hzn_master[i])) != 0) & (
    is.na(dat$hzn_top[i])==F)) {
    if ((length(which(horizonhold==dat$hzn_desgn[i])) != 0) &
          (dat$mid[i]<=mindepth)) {
      hhold <- c(hhold,i)
    }
  }
}
noA <- length(hhold)
d <- dat[hhold,]
names(d)[which(names(d)=="sand_tot_psa")] <- "sand"
names(d)[which(names(d)=="silt_tot_psa")] <- "silt"
names(d)[which(names(d)=="clay_tot_psa")] <- "clay"
names(d)[which(names(d)=="tex_psda")] <- "txt"
txt <- tolower(as.character(d$txt))
ltxt <- levels(as.factor(txt))
for(i in 1:length(ltxt)) {
  if (ltxt[i]=="cos") txt[which(txt=="cos")] <- "s"
  if (ltxt[i]=="cosl") txt[which(txt=="cosl")] <- "sl"
  if (ltxt[i]=="fs") txt[which(txt=="fs")] <- "s"
  if (ltxt[i]=="fsl") txt[which(txt=="fsl")] <- "sl"
  if (ltxt[i]=="lcos") txt[which(txt=="lcos")] <- "ls"
  if (ltxt[i]=="lfs") txt[which(txt=="lfs")] <- "ls"
  if (ltxt[i]=="lvfs") txt[which(txt=="lvfs")] <- "ls"
  if (ltxt[i]=="vfsl") txt[which(txt=="vfsl")] <- "sl"
}
d$txt <- as.factor(txt)
rm(list=c("txt","ltxt"))
#
d <- d[-which(is.na(d$oc)),] # removes SOC NA values
#
dA <- d # saves the data frame for this selection
#
# Save the lowest point in Florida for removal below (i.e., FLremove==T
# is set above)
FLmin <- min(d$lat[which(d$state_admindiv_code=="FL")])
#
#
# Interpolated maps
# _________________
rerun <- T # set to false if the climate plots are not needed
if (rerun==T) {
  # Climate
  pdf(file="ClimateMaps.pdf",height=6,width=3)
  par(mar=c(0,0,0,0),mfrow=c(3,1),oma=c(2,2,0,0))
  setwd(sname)
  source('climateinterpolation.R')
  setwd(fname)
  dev.off()
}
# Porosities
# Set the scale after running both the A and the B horizons
minvalueep <- 0.09; maxvalueep <- 0.3
minvaluewhc <- 0.06; maxvaluewhc <- 0.22
pdf(file=paste(mtlt,"InterplotatedMaps.pdf",sep="_"),height=6,width=6)
par(mar=c(0,0,0,0),mfrow=c(3,2),oma=c(2,2,0,0))
setwd(sname)
source('interpolation.R')
setwd(fname)
dev.off()
#
#
# Calculate aggregated data for the climate dataframe
# ---------------------------------------------------
#
#  Finds and calculates the mean EP and residual EP for the closest station
rerun <- T # set to false if ep points to the nearest station are already 
#   found
if (rerun==T) {
  # Checks to see if the climate data.frame already contains these
  # variables and if so, removes them
  climatehold <- which((names(climate)=="Aep") |
                         ((names(climate)=="Aresep") |
                            ((names(climate)=="Atp") |
                               ((names(climate)=="Afc") |
                                  ((names(climate)=="Awp") |
                                     ((names(climate)=="Asand") |
                                        ((names(climate)=="Aclay")) | 
                                        ((names(climate)=="Arestp")) | # NEW
                                        ((names(climate)=="Aresfc")) | # NEW
                                        (names(climate)=="Areswp"))))))) # NEW
  if (length(climatehold)!=0) {
    climate <- climate[,-climatehold] # removes the offending columns
  }
  #
  # Load the package 'fields' (See note in the preamble above)
  require(fields)
  #
  # Create a distance matrix to find the minimum distance for the 
  #   A horizons
  aquichold <- which(dA$moistreg=="AQUIC") # for removing aquic moisture 
  # regime samples
  dAhold <- dA # saves the dA data with the aquic moisture regime
  dA <- dA[-aquichold,] # removes the aquic moisture regime samples
  cl <- as.matrix(data.frame(lon=climate$long,lat=climate$lat))
  da <- as.matrix(data.frame(lon=dA$long,lat=dA$lat))
  distmat <- rdist.earth(cl,da) # rows are the row index for cl; cols are the
  #   row index for da
  #
  # Run a loop to record the nearest station to each climate station
  stmat <- matrix(NA,nrow=nrow(distmat),ncol=ncol(distmat))
  for (j in 1:ncol(distmat)) {
    stmat[which.min(distmat[,j]),j] <- 1
  }
  #
  # Run a loop to calculate mean effective porosity for the A horizons
  steff <- NULL
  for (i in 1:nrow(stmat)) {
    hold <- which(stmat[i,]==1)
    if (length(hold)>0) {
      steff <- c(steff,mean(dA$ep[hold]))
    } else {
      steff <- c(steff,NA)
    }
  }
  climate <- data.frame(climate,Aep=steff)
  #
  # Run a loop to calculate mean total porosity for the A horizons
  steff <- NULL
  for (i in 1:nrow(stmat)) {
    hold <- which(stmat[i,]==1)
    if (length(hold)>0) {
      steff <- c(steff,mean(dA$tp[hold]))
    } else {
      steff <- c(steff,NA)
    }
  }
  climate <- data.frame(climate,Atp=steff)
  #
  # Run a loop to calculate mean field capacity for the A horizons
  steff <- NULL
  for (i in 1:nrow(stmat)) {
    hold <- which(stmat[i,]==1)
    if (length(hold)>0) {
      steff <- c(steff,mean(dA$fc[hold]))
    } else {
      steff <- c(steff,NA)
    }
  }
  climate <- data.frame(climate,Afc=steff)
  #
  # Run a loop to calculate mean wilting point for the A horizons
  steff <- NULL
  for (i in 1:nrow(stmat)) {
    hold <- which(stmat[i,]==1)
    if (length(hold)>0) {
      steff <- c(steff,mean(dA$wp[hold]))
    } else {
      steff <- c(steff,NA)
    }
  }
  climate <- data.frame(climate,Awp=steff)
  #
  # Run a loop to calculate mean sand content for the A horizons
  steff <- NULL
  for (i in 1:nrow(stmat)) {
    hold <- which(stmat[i,]==1)
    if (length(hold)>0) {
      steff <- c(steff,mean(dA$sand[hold]))
    } else {
      steff <- c(steff,NA)
    }
  }
  climate <- data.frame(climate,Asand=steff)
  #
  # Run a loop to calculate mean clay content for the A horizons
  steff <- NULL
  for (i in 1:nrow(stmat)) {
    hold <- which(stmat[i,]==1)
    if (length(hold)>0) {
      steff <- c(steff,mean(dA$clay[hold]))
    } else {
      steff <- c(steff,NA)
    }
  }
  climate <- data.frame(climate,Aclay=steff)
}
# Save the dA data without the aquic moisture regime
dAn <- dA
# Put the aquic moisture regime samples back into the dA data.frame
dA <- dAhold
rm(dAhold)
#
# ---------------------------------------------------------------------------
