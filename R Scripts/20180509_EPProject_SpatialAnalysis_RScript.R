# ---------------------------------------------------------------------------
# Title: R log for effective porosity manuscript
# Author: D.R. Hirmas
# Explanation: Subset data to provide to Ruth Kerry for spatial data
# analysis in response to the second Nature review
# ---------------------------------------------------------------------------
#
#
# ---------------------------------------------------------------------------
# Set the working directory, load packages, and read in data
# ---------------------------------------------------------------------------
# Save the working directory and subfolder directories:
setwd("..") # assumes the RStudio file is in a separate folder
wname <- getwd() # parent folder for the project
#
# Save the working directory where the R Script files are stored:
sname <- paste(wname,"R Scripts",sep="/")
#
# Save the working directory where the data files are stored:
dname <- paste(wname,"Data",sep="/")
#
# Change the working directory where the output files will be stored:
fname <- paste(wname,"Plots",sep="/")
oname <- paste(wname,"Output Files",sep="/")
#
# Save the working directory where the function files are stored:
funname <- paste(wname,"Functions",sep="/")
#
# Retrieve the latest function file name (assumes all functions are stored
# in one file, that 'Functions' is somewhere in the file name exactly, and
# that the file name begins with a date in the following format:
# YearMonthDay; this formate comprises a total of 8 characters)
setwd(funname)
funhold <- list.files(pattern='Functions')
funhold <- funhold[which.max(as.numeric(substr(funhold,1,8)))]
setwd(wname)
#
# Read in the data files:
setwd(dname)
dat <- read.delim("nrcsdata.txt",header=T) # holds the original soil
# property information; Note: this data.frame is rewritten in the 
# "Find the nearest station and apply those timing and magnitude values" 
# section below if rerun is set to T
nrcskoep <- read.delim("NRCS2.txt",header=T) # NRCS2 holds the Koeppen, 
# moisture regime, and temperature regime data; file created by Rob 
# Miskewitz from the shape files obtained from Cassie Wilson
#
# Weather station (US Historical Climatology Network)
stn <- read.delim("stn_sm_full.txt",header=T) # holds the USHCN weather
# station data converted to timing and magnitude parameters following
# D’Odorico et al. (2000)
ppt <- read.table("1951-2011_ave_ppt.txt",sep=" ") # holds the mean annual
# precipitation record for each station (mm)
names(ppt) <- c("stnid","lat","long","ppt")
temp <- read.table("1951-2011_ave_temps.txt",sep=" ",skip=1) # holds the 
# mean max and min temperatures for each station (C)
names(temp) <- c("stnid","lat","long","maxtemp","mintemp")
freeze <- read.table("1951-2011_freeze.txt",sep=" ",skip=1) # holds the 
# number of times that the min or max temp dropped below freezing for each
# year of the record
names(freeze) <- c("stnid","lat","long","year","mintempfreezedays",
                   "mintempfreezedays")
#
setwd(sname) # need to run the following scripts to access the data
dnamehold <- dname # change of variables to accomodate the use of dname
# in the following scripts
# Climate interpolation (PRISM)
source('prismdatareadin.R') # reads in the PRISM data for precipitation
# temperature, dewpoint temperature, min and max vapor pressure deficit,
# min and max temperature, and elevation (4 km grid) interpolated from
# 1981-2010 normals. Names of the new data.frames are: prismdewtemp (C), 
# prismmaxtemp (C), preismmintemp (C), prismminvpd (hPa), prismmaxvpd (hPa),
# prismppt (mm), prismtemp (C), and elev (m)
setwd(sname) # resets the working directory as it is changed in the
# previous source code
dname <- dnamehold # resets the dname as it is changed in the previous
# source code
#
# Climate prediction (Climate Wizard)
source('climatewizarddatareadin.R') # reads in the Climate Wizard data for 
# end of century (2080) precipitation under the medium A1B senario from an 
# average of climate model predictions.
setwd(sname) # resets the working directory as it is changed in the
# previous source code
dname <- dnamehold # resets the dname as it is changed in the previous
# source code
#
# Climate prediction (CESM Data from Nate)
source('cesmdatareadin.R') # reads in CESM data for end of century (2100)
# under a middle-of-the road scenario
dname <- dnamehold # resets the dname as it is changed in the previous
# source code
setwd(wname) # resets the working directory as it is changed in the
# previous source code
#
# Load the following packages:
library(maps); library(gstat); library(MASS); library(raster); 
library(matlab); library(segmented); library(MVN); library(lmodel2)
library(car); library(lattice); library(spatial); require(colorRamps);
library(RColorBrewer); library(plyr); library(grid); require(normwhn.test)
require(geosphere) # needed for the distHaversine function in the new code
# require(fields) # this will mask matlab so it is not added until needed 
# below
#
# Load the following predefined functions:
setwd(funname)
source(funhold)
setwd(fname)
#
# Change to pdf or eps for the final figures:
# [Use in the wrapper:]
# if (filetype==T) {
#   pdf(file=paste(gsub("-", "", as.character(Sys.Date())),
#     '_FileName.pdf',sep=''), 
#     width=PageWidthinInches, 
#     height=PageHeightinInches)
# } else {
#   eps(file=paste(gsub("-", "", as.character(Sys.Date())),
#      '_FileName.eps',sep=''), 
#      width=PageWidthinInches, 
#      height=PageHeightinInches)
# }
figuretype <- "pdf"
# 
# [Use to output RData files:]
# setwd(oname)
# save(list='RObject(s)',file=paste(gsub("-", "", as.character(Sys.Date())),
#    '_FileName.RData',sep=""))
# setwd(wname)
# [or]
# setwd(oname)
# write.table(dAn,file=paste(gsub("-", "", as.character(Sys.Date())),
#             '_FileName.csv',sep=""),sep=",",row.names=F)
# setwd(wname)
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Convert the lat long data to decimal degrees and clean up the data.frame
# ---------------------------------------------------------------------------
# Convert the lat long data
lat <- dat$latitude_degrees + dat$latitude_minutes/60 + 
  dat$latitude_seconds/3600
long <- dat$longitude_degrees + dat$longitude_minutes/60 + 
  dat$longitude_seconds/3600
#
# Change the column headings and add the lat long data
dat <- dat[,-which(names(dat)=="latitude_minutes")]
dat <- dat[,-which(names(dat)=="latitude_seconds")]
dat <- dat[,-which(names(dat)=="longitude_minutes")]
dat <- dat[,-which(names(dat)=="longitude_seconds")]
names(dat)[which(names(dat)=="latitude_degrees")] <- "lat"
names(dat)[which(names(dat)=="longitude_degrees")] <- "long"
dat$lat <- lat; dat$long <- -long
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Add effective and total porosity and water holding capacity
# ---------------------------------------------------------------------------
# Calculate total porosity
tp <- (1-dat$db_13b/2.65)
#
# Calculate field capacity water content
fc <- dat$w3cld/100*dat$db_13b # converts to volumetric
#
# Calculate effective porosity
ep <- tp-fc
#
# Calculate wilting point water content
wp <- dat$w15l2/100*dat$db_od # converts to volumetric
#
# Calculate available-water holding capacity
#whc <- fc-dat$w15l2/100
whc <- fc-wp
#
# Add these to the data.frame
dat <- data.frame(dat,tp,fc,wp,ep,whc)
#
# Remove those data points where wp is higher than fc
phold <- which(dat$whc>0)
dat <- dat[phold,]
#
# Remove those data points outside the mainland
phold <- which(((dat$lat < 50) & (dat$lat > 24)) & 
                 ((dat$long < -66) & (dat$long > -125)))
dat <- dat[phold,]
# 
# Remove outlier data (trims bottom 1% and top 99% for 3 and 15 bar water
# contents, total porosity, and effective porosity and removes nonsensical
# data such as samples with 15 bar water contents larger than the 3 bar or
# 3 bar water contents larger than total porosity)
lw3 <- quantile(dat$fc,probs=0.01); uw3 <- quantile(dat$fc,probs=0.99)
lw15 <- quantile(dat$wp,probs=0.01); uw15 <- quantile(dat$wp,probs=0.99)
ltp <- quantile(dat$tp,probs=0.01); utp <- quantile(dat$tp,probs=0.99)
lep <- quantile(dat$ep,probs=0.01); uep <- quantile(dat$ep,probs=0.99)
lwhc <- quantile(dat$whc,probs=0.01); uwhc <- quantile(dat$whc,probs=0.99)
whold <- which((((((dat$fc < uw3) & (dat$fc > lw3)) & 
                    ((dat$wp < uw15) & (dat$wp > lw15))) & 
                   (((dat$tp < utp) & (dat$tp > ltp)) & 
                      ((dat$ep < uep) & (dat$ep > lep)))) & ((dat$tp > dat$fc) & 
                                                               (dat$fc > dat$wp))) & ((dat$whc > lwhc) & (dat$whc < uwhc)))
dat <- dat[whold,]
#
# Remove the w3cld and w15l2 columns since those are gravimetric
phold <- which((colnames(dat)=='w3cld') | (colnames(dat)=='w15l2'))
dat <- dat[,-phold]
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Weather station data
# ---------------------------------------------------------------------------
# Create a new dataframe with lat, long, mag, and timing
alpha <- stn[,"POISSON"]
lambda <- stn[,"EXPONENTIAL"]
id <- stn[,"STNID"]
al <- data.frame(stnid=id,lat=stn[,"LATITUDE"],long=stn[,"LONGITUDE"],
                 magnitude=1/lambda,timing=alpha) # this is correct
#
# Merge all the data together
climate <- merge(temp,ppt,by.x="stnid",by.y="stnid")
climate <- merge(climate,al,by.x="stnid",by.y="stnid")
hold <- which((names(climate)=="lat.y") | (names(climate)=="long.y") |
                (names(climate)=="lat") | (names(climate)=="long"))
climate <- climate[,-hold]
names(climate)[which(names(climate)=="lat.x")] <- "lat"
names(climate)[which(names(climate)=="long.x")] <- "long"
#
# Summarize the freeze data for each station and add to climate
freezehold <- aggregate(cbind(maxtempfreezedays,mintempfreezedays)~
                          stnid,freeze,mean)
climate <- data.frame(climate,maxtempfreezedays=
                        freezehold[,"maxtempfreezedays"],
                      mintempfreezedays=
                        freezehold[,"mintempfreezedays"])
rm(freezehold)
#
#
# Add the PRISM data to the climate data.frame
climate <- data.frame(climate,prismppt=rep(NA,nrow(climate)),
                      prismtemp=rep(NA,nrow(climate)), 
                      prismmintemp=rep(NA,nrow(climate)),
                      prismmaxtemp=rep(NA,nrow(climate)),
                      prismdewtemp=rep(NA,nrow(climate)),
                      #prismminvpd=rep(NA,nrow(climate)),
                      prismjjavpd=rep(NA,nrow(climate)), # mean JJA daily max vpd 
                      prismmaxvpd=rep(NA,nrow(climate)),
                      elev=rep(NA,nrow(climate)))
for (i in 1:nrow(climate)) {
  cr <- max(which(colref < climate[i,"long"])) # colref comes from 
  # the prismdatareadin.R file
  nr <- min(which(rowref < climate[i,"lat"])) # same for rowref
  climate[i,11:18] <- c(prismppt[nr,cr],prismtemp[nr,cr],prismmintemp[nr,cr],
                        prismmaxtemp[nr,cr],prismdewtemp[nr,cr],
                        #prismminvpd[nr,cr],
                        prismjjavpd[nr,cr],prismmaxvpd[nr,cr],elev[nr,cr])
  # Note: Careful! Hardcoded the column placement of the prism data
} # Check this bracket...this may be an error!!
#
# Add the climate wizard data to the climate data.frame
climate <- data.frame(climate,wizardppt=rep(NA,nrow(climate)))
for (i in 1:nrow(climate)) {
  cr <- max(which(wcolref < climate[i,"long"])) # wcolref comes from 
  # the climatewizarddatareadin.R file
  nr <- min(which(wrowref < climate[i,"lat"])) # same for wrowref
  if(cr == length(wcolref)) cr <- length(wcolref)-1
  if(nr == length(wrowref)) nr <- length(wrowref)-1
  climate[i,ncol(climate)] <- wizardppt[nr,cr]
}
climate$wizardppt <- climate$wizardppt*25.4 # converts from inches to mm
#
# Add the CESM climate data to the climate data.frame
climate <- data.frame(climate,cesmppt=rep(NA,nrow(climate)))
for (i in 1:nrow(climate)) {
  cr <- max(which(ccolref < climate[i,"long"])) # ccolref comes from 
  # the cesmdatareadin.R file
  nr <- min(which(crowref < climate[i,"lat"])) # same for crowref
  if(cr == length(ccolref)) cr <- length(ccolref)-1
  if(nr == length(crowref)) nr <- length(crowref)-1
  climate[i,ncol(climate)] <- cesmppt[nr,cr]
}
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Find the nearest station and apply those climate values
# ---------------------------------------------------------------------------
rerun <- F # set to false if nearest station is already found
if (rerun==T) {
  #
  # Load the package 'fields'
  require(fields)
  #
  # Create a distance matrix to find the minimum distance
  cl <- as.matrix(data.frame(lon=climate$long,lat=climate$lat))
  da <- as.matrix(data.frame(lon=dat$long,lat=dat$lat))
  distmat <- rdist.earth(da,cl) # rows are the row index for da; cols are the
  # 	row index for cl
  stnord <- NULL
  for (i in 1:nrow(distmat)) stnord <- c(stnord,which.min(distmat[i,]))
  stnid <- climate$stnid[stnord] # added 2/1/18
  mag <- climate$magnitude[stnord]
  timing <- climate$timing[stnord]
  ppt <- climate$ppt[stnord]
  mintemp <- climate$mintemp[stnord]
  maxtemp <- climate$maxtemp[stnord]
  maxtempfreezedays <- climate$maxtempfreezedays[stnord]
  mintempfreezedays <- climate$mintempfreezedays[stnord]
  dat <- data.frame(dat,stnid,mag,timing,ppt,mintemp,maxtemp,
                    maxtempfreezedays,mintempfreezedays)
  #
  # Separately add the PRISM data to each profile
  startcolnum <- ncol(dat)+1 # the starting column to put the data below
  endcolnum <- startcolnum+7 # the ending column to put the data below
  dat <- data.frame(dat,prismppt=rep(NA,nrow(dat)),
                    prismtemp=rep(NA,nrow(dat)),
                    prismmintemp=rep(NA,nrow(dat)),
                    prismmaxtemp=rep(NA,nrow(dat)),
                    prismdewtemp=rep(NA,nrow(dat)),
                    prismjjavpd=rep(NA,nrow(dat)),
                    prismmaxvpd=rep(NA,nrow(dat)),
                    elev=rep(NA,nrow(dat)))
  for (i in 1:nrow(dat)) {
    cr <- max(which(colref < dat[i,"long"])) # colref comes from 
    # the prismdatareadin.R file
    nr <- min(which(rowref < dat[i,"lat"])) # same for rowref
    dat[i,startcolnum:endcolnum] <- c(prismppt[nr,cr],
                                      prismtemp[nr,cr],
                                      prismmintemp[nr,cr],
                                      prismmaxtemp[nr,cr],
                                      prismdewtemp[nr,cr],
                                      prismjjavpd[nr,cr],
                                      prismmaxvpd[nr,cr],
                                      elev[nr,cr])
  }
  corcolnum <- ncol(dat)+1 # the correct column to put the data below
  dat <- data.frame(dat,cesmppt=rep(NA,nrow(dat)))
  for (i in 1:nrow(dat)) {
    cr <- max(which(ccolref < dat[i,"long"])) # ccolref comes from 
    # the cesmdatareadin.R file
    nr <- min(which(crowref < dat[i,"lat"])) # same for crowref
    dat[i,corcolnum] <- cesmppt[nr,cr]
  }
  #
  # Write the file to read in on another machine
  setwd(dname) # saves it in the Data folder
  write.table(dat,file="datcl2.txt") # modified 2/1/18
  setwd(fname)
} else {setwd(dname); dat <- read.table("datcl2.txt"); setwd(fname)}
# Note: This else statement reads in the new datcl2.txt data (includes
# stnid information)
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Utilize the Koeppen, moisture regime, and the temperature regime data
# ---------------------------------------------------------------------------
# Obtain the Koeppen, moisture regime, and temp regime classification from
# nrcskoep data.frame read in earlier
koeppen <- NULL
moistreg <- NULL
tempreg <- NULL
for (i in 1:nrow(dat)) {
  hold <- which(((nrcskoep$pedon_key==dat$pedon_key[i]) & 
                   (nrcskoep$layer_key==dat$layer_key[i])) & 
                  (nrcskoep$site_key==dat$site_key[i]))
  if (length(hold)==1) {
    koeppen <- c(koeppen,as.character(nrcskoep[hold,"Koeppen__1"]))
    moistreg <- c(moistreg,as.character(nrcskoep[hold,"REGIME_1"]))
    tempreg <- c(tempreg,as.character(nrcskoep[hold,"REGIME"]))
  } else {koeppen <- c(koeppen,NA); moistreg <- c(moistreg,NA);
  tempreg <- c(tempreg,NA)}
}
#
# Add this data to the dat data.frame
dat <- data.frame(dat,koeppen,moistreg,tempreg)
rm(list=c('nrcskoep','koeppen','moistreg','tempreg')) # removes uneeded
# files
#
# Add the midpoint depth
mid <- NULL
for (i in 1:nrow(dat)) {
  mid <- c(mid,mean(c(dat$hzn_top[i],dat$hzn_bot[i])))
}
dat <- data.frame(dat,mid)
#
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Examine the relationship between climate and porosity for A horizons (see
# the selection of A horizons below) that have their mid horizon depth  
# within the upper 25 cm
# ---------------------------------------------------------------------------
# Set the following to F if the Florida interpolation removal (in the case
# of no points) is not needed
FLremove <- T
#
setwd(sname)
source('Ahorizon.R') # this generates quite a number of plots
# Note: The lattice plots do not work when calling them remotely in a 
# source file. To get those to work, they must be run directly from the
# source file. These plots are the EP or residual EP against the climate
# variables for each of the 12 soil textural classes.
#
# Summarize the results of the regressions
Apptlm2 # residual effective porostiy vs MAP
summary(Apptlm)
Amaglm2 # residual effective porostiy vs event magnitude
summary(Amaglm)
Avpdlm2 # residual effective porostiy vs vapor pressure deficit
summary(Avpdlm)
Atemplm2 # residual effective porostiy vs MAT
summary(Atemplm)
Afreezelm2 # residual effective porostiy vs freeze frequency
summary(Afreezelm)
Aeemtlm2 # residual effective porostiy vs EEMT
summary(Aeemtlm)
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Examine the relationship between climate and porosity for Ap horizons (see
# the selection of Ap horizons below) that have their mid horizon depth  
# within the upper 25 cm
# ---------------------------------------------------------------------------
# Set the following to F if the Florida interpolation removal (in the case
# of no points) is not needed
FLremove <- T
#
setwd(sname)
source('Aphorizon.R') # this generates quite a number of plots
# Note: The lattice plots do not work when calling them remotely in a 
# source file. To get those to work, they must be run directly from the
# source file. These plots are the EP or residual EP against the climate
# variables for each of the 12 soil textural classes.
#
# Summarize the results of the regressions
Appptlm2 # residual effective porostiy vs MAP
summary(Appptlm)
Apmaglm2 # residual effective porostiy vs event magnitude
summary(Apmaglm)
Apvpdlm2 # residual effective porostiy vs vapor pressure deficit
summary(Apvpdlm)
Aptemplm2 # residual effective porostiy vs MAT
summary(Aptemplm)
Apfreezelm2 # residual effective porostiy vs freeze frequency
summary(Apfreezelm)
Apeemtlm2 # residual effective porostiy vs EEMT
summary(Apeemtlm)
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Examine the relationship between climate and porosity for B horizons (see
# the selection of B horizons below) that have their mid horizon depth 
# between 25 and 100 cm
# ---------------------------------------------------------------------------
# Set the following to F if the Florida interpolation removal (in the case
# of no points) is not needed
FLremove <- T
#
setwd(sname)
source('Bhorizon.R') # this generates quite a number of plots
# Note: The lattice plots do not work when calling them remotely in a 
# source file. To get those to work, they must be run directly from the
# source file. These plots are the EP or residual EP against the climate
# variables for each of the 12 soil textural classes.
#
# Summarize the results of the regressions
Bpptlm2 # residual effective porostiy vs MAP
summary(Bpptlm)
Bmaglm2 # residual effective porostiy vs event magnitude
summary(Bmaglm)
Bvpdlm2 # residual effective porostiy vs vapor pressure deficit
summary(Bvpdlm)
Btemplm2 # residual effective porostiy vs MAT
summary(Btemplm)
Bfreezelm2 # residual effective porostiy vs freeze frequency
summary(Bfreezelm)
Beemtlm2 # residual effective porostiy vs EEMT
summary(Beemtlm)
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Add climate station lat long to each of the dataframes
# ---------------------------------------------------------------------------
climateholddf <- climate[,c("stnid","lat","long")]
names(climateholddf) <- c("stnid","stnlat","stnlong")
dAn <- merge(dAn,climateholddf,by='stnid',all.x=T)
dApn <- merge(dApn,climateholddf,by='stnid',all.x=T)
dBn <- merge(dBn,climateholddf,by='stnid',all.x=T)
rm(list='climateholddf')
#
# Output the dataframes as separate files
setwd(oname)
write.table(dAn,file=paste(gsub("-", "", as.character(Sys.Date())),
                           '_EPProject_AHorizons.csv',sep=""),
            sep=",",row.names=F)
write.table(dApn,file=paste(gsub("-", "", as.character(Sys.Date())),
                           '_EPProject_ApHorizons.csv',sep=""),
            sep=",",row.names=F)
write.table(dBn,file=paste(gsub("-", "", as.character(Sys.Date())),
                           '_EPProject_BHorizons.csv',sep=""),
            sep=",",row.names=F)
setwd(fname)
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Create the interpolated resepx maps that incorporate the spatial model
# in the calculation of residual effective porosity; these resepx values
# come from Ruth Kerry (added 3/22/18)
# ---------------------------------------------------------------------------
# Set the following to F if the Florida interpolation removal (in the case
# of no points) is not needed
FLremove <- T
#
require(fields)
filetype <- "pdf"
PageWidthinInches <- 6
PageHeightinInches <- 6
setwd(dname)
d <- read.csv('20180314_EPProject_AHorizons_RESEPCorrectedwithaSpatialModel_aRK.csv')
dAx <- d
setwd(fname)
mtlt <- "AHorizonswithMidpointDepthswithin25cm"
Filename <- paste(mtlt,"ResidualInterplotatedMaps",
                  "SpatialModelResidualEP",sep="_")
minvalueresep <- -0.08; maxvalueresep <- 0.08 # use for a divergent ramp
minvaluereswhc <- -0.05; maxvaluereswhc <- 0.04
if (filetype=="pdf") {
  pdf(file=paste(gsub("-", "", as.character(Sys.Date())),
                 "_",Filename,'.pdf',sep=''), 
      width=PageWidthinInches, 
      height=PageHeightinInches)
} else { # .eps
  eps(file=paste(gsub("-", "", as.character(Sys.Date())),
                 "_",Filename,'.eps',sep=''), 
      width=PageWidthinInches,
      height=PageHeightinInches)
}
par(mar=c(0,0,0,0),mfrow=c(3,2),oma=c(2,2,0,0))
setwd(sname)
source('resxinterpolation.R')
dev.off()
setwd(dname)
d <- read.csv('20180314_EPProject_ApHorizons_RESEPCorrectedwithaSpatialModel_aRK.csv')
dApx <- d
setwd(fname)
mtlt <- "ApHorizonswithMidpointDepthswithin25cm"
Filename <- paste(mtlt,"ResidualInterplotatedMaps",
                  "SpatialModelResidualEP",sep="_")
if (filetype=="pdf") {
  pdf(file=paste(gsub("-", "", as.character(Sys.Date())),
                 "_",Filename,'.pdf',sep=''), 
      width=PageWidthinInches, 
      height=PageHeightinInches)
} else { # .eps
  eps(file=paste(gsub("-", "", as.character(Sys.Date())),
                 "_",Filename,'.eps',sep=''), 
      width=PageWidthinInches,
      height=PageHeightinInches)
}
par(mar=c(0,0,0,0),mfrow=c(3,2),oma=c(2,2,0,0))
setwd(sname)
source('resxinterpolation.R')
dev.off()
setwd(dname)
d <- read.csv('20180314_EPProject_BHorizons_RESEPCorrectedwithaSpatialModel_aRK.csv')
dBx <- d
setwd(fname)
mtlt <- "BHorizonswithMidpointDepthswithin25cm"
Filename <- paste(mtlt,"ResidualInterplotatedMaps",
                  "SpatialModelResidualEP",sep="_")
if (filetype=="pdf") {
  pdf(file=paste(gsub("-", "", as.character(Sys.Date())),
                 "_",Filename,'.pdf',sep=''), 
      width=PageWidthinInches, 
      height=PageHeightinInches)
} else { # .eps
  eps(file=paste(gsub("-", "", as.character(Sys.Date())),
                 "_",Filename,'.eps',sep=''), 
      width=PageWidthinInches,
      height=PageHeightinInches)
}
par(mar=c(0,0,0,0),mfrow=c(3,2),oma=c(2,2,0,0))
setwd(sname)
source('resxinterpolation.R')
dev.off()
setwd(fname)
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Aggregate the RESEPx data to the weather stations in the same way that that 
# was done for the original Nature submission and send to Ruth to run spatial 
# regression (3/23/18); modified again (5/9/18)
# ---------------------------------------------------------------------------
# Aggregate the RESEPx data and added it to the climate data.frame
setwd(sname)
source('Axhorizon.R')
source('Apxhorizon.R')
source('Bxhorizon.R')
#
# Remove outliers using an adjusted Mahalanobis distance method following 
# Korkmaz et al. (2014) (this is the way that the original Nature submission
# had removed the outliers):
#
# Function to identify bivariate outliers using the adjusted Mahalanobois 
# distance procedure
outlierID <- function(df) {
  # Creates a logical vector to test to see if the point is a bivariate
  # outlier; vector corresponds to the order of the rows in the input df
  # df = data.frame; only the first two columns will be used to detect the
  # outliers with; those first two columns need to be numeric
  #
  # Please note: there is a bizzare problem with mvOutlier where if the 
  # row.names of the dataframe are NULL from 1 to nrow(df), the algorithm
  # will designated the number of outliers correctly but will remove the
  # wrong ones in the resulting dfout$newData dataframe. Solution is to 
  row.names(df) <- 1:nrow(df)
  dfout <- mvOutlier(df,qqplot=F,method="adj.quan")
  dfo <- dfout$outlier # data.frame indicating the outliers
  dfnd <- dfout$newData # df without the outliers
  ifcheck <- length(which(dfo$Outlier==TRUE)) > 0
  if (ifcheck) {
    out <- rep(FALSE,nrow(df))
    dfindex <- as.numeric(row.names(dfo[which(dfo$Outlier==TRUE),]))
    out[dfindex] <- TRUE
  } else {
    out <- rep(FALSE,nrow(df))
  }
  out # TRUE indicates outlier; FALSE indicates not an outlier
}
#
#
# Removes outliers for the resepx vs ppt columns
dhold <- data.frame(ppt=climate$ppt,resepx=climate$Aresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Aresepxpptout <- outlierID(dhold)
Aresepxpptout <- data.frame(outlier=outlierID(dhold),
                            stnid=dhold$stnid)
Aresepxpptout <- base::merge(climate,
                             Aresepxpptout,
                             by="stnid",
                             all.x=TRUE)$outlier
#
dhold <- data.frame(ppt=climate$ppt,resepx=climate$Apresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Apresepxpptout <- outlierID(dhold)
Apresepxpptout <- data.frame(outlier=outlierID(dhold),
                             stnid=dhold$stnid)
Apresepxpptout <- base::merge(climate,
                              Apresepxpptout,
                              by="stnid",
                              all.x=TRUE)$outlier
#
dhold <- data.frame(ppt=climate$ppt,resepx=climate$Bresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Bresepxpptout <- data.frame(outlier=outlierID(dhold),
                            stnid=dhold$stnid)
Bresepxpptout <- base::merge(climate,
                             Bresepxpptout,
                             by="stnid",
                             all.x=TRUE)$outlier
#
# Removes outliers for the resepx vs prismjjavpd columns
dhold <- data.frame(prismjjavpd=climate$prismjjavpd,
                    resepx=climate$Aresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Aresepxprismjjavpdout <- data.frame(outlier=outlierID(dhold),
                                    stnid=dhold$stnid)
Aresepxprismjjavpdout <- base::merge(climate,
                                     Aresepxprismjjavpdout,
                                     by="stnid",
                                     all.x=TRUE)$outlier
#
dhold <- data.frame(prismjjavpd=climate$prismjjavpd,
                    resepx=climate$Apresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Apresepxprismjjavpdout <- data.frame(outlier=outlierID(dhold),
                                     stnid=dhold$stnid)
Apresepxprismjjavpdout <- base::merge(climate,
                                      Apresepxprismjjavpdout,
                                      by="stnid",
                                      all.x=TRUE)$outlier
#
dhold <- data.frame(prismjjavpd=climate$prismjjavpd,
                    resepx=climate$Bresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Bresepxprismjjavpdout <- data.frame(outlier=outlierID(dhold),
                                    stnid=dhold$stnid)
Bresepxprismjjavpdout <- base::merge(climate,
                                     Bresepxprismjjavpdout,
                                     by="stnid",
                                     all.x=TRUE)$outlier
#
# Removes outliers for the resepx vs magnitude columns
dhold <- data.frame(magnitude=climate$magnitude,resepx=climate$Aresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Aresepxmagnitudeout <- data.frame(outlier=outlierID(dhold),
                                  stnid=dhold$stnid)
Aresepxmagnitudeout <- base::merge(climate,
                                   Aresepxmagnitudeout,
                                   by="stnid",
                                   all.x=TRUE)$outlier
#
dhold <- data.frame(magnitude=climate$magnitude,resepx=climate$Apresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Apresepxmagnitudeout <- data.frame(outlier=outlierID(dhold),
                                   stnid=dhold$stnid)
Apresepxmagnitudeout <- base::merge(climate,
                                    Apresepxmagnitudeout,
                                    by="stnid",
                                    all.x=TRUE)$outlier
#
dhold <- data.frame(magnitude=climate$magnitude,resepx=climate$Bresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Bresepxmagnitudeout <- data.frame(outlier=outlierID(dhold),
                                  stnid=dhold$stnid)
Bresepxmagnitudeout <- base::merge(climate,
                                   Bresepxmagnitudeout,
                                   by="stnid",
                                   all.x=TRUE)$outlier
#
# Removes outliers for the resepx vs mintempfreezedays columns
dhold <- data.frame(mintempfreezedays=climate$mintempfreezedays,
                    resepx=climate$Aresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Aresepxmintempfreezedaysout <- data.frame(outlier=outlierID(dhold),
                                          stnid=dhold$stnid)
Aresepxmintempfreezedaysout <- base::merge(climate,
                                           Aresepxmintempfreezedaysout,
                                           by="stnid",
                                           all.x=TRUE)$outlier
#
dhold <- data.frame(mintempfreezedays=climate$mintempfreezedays,
                    resepx=climate$Apresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Apresepxmintempfreezedaysout <- data.frame(outlier=outlierID(dhold),
                                           stnid=dhold$stnid)
Apresepxmintempfreezedaysout <- base::merge(climate,
                                            Apresepxmintempfreezedaysout,
                                            by="stnid",
                                            all.x=TRUE)$outlier
#
dhold <- data.frame(mintempfreezedays=climate$mintempfreezedays,
                    resepx=climate$Bresepx,
                    stnid=climate$stnid)
dhold <- dhold[complete.cases(dhold),] # removes rows with NAs
Bresepxmintempfreezedaysout <- data.frame(outlier=outlierID(dhold),
                                          stnid=dhold$stnid)
Bresepxmintempfreezedaysout <- base::merge(climate,
                                           Bresepxmintempfreezedaysout,
                                           by="stnid",
                                           all.x=TRUE)$outlier
#
climatehold <- which((names(climate)=="Aresepxpptout") |
                       ((names(climate)=="Apresepxpptout") |
                       ((names(climate)=="Bresepxpptout") |
                       ((names(climate)=="Aresepxprismjjavpdout") |
                       ((names(climate)=="Apresepxprismjjavpdout") |
                       ((names(climate)=="Bresepxprismjjavpdout") |
                       ((names(climate)=="Aresepxmagnitudeout") |
                       ((names(climate)=="Apresepxmagnitudeout") |
                       ((names(climate)=="Bresepxmagnitudeout") |
                       ((names(climate)=="Aresepxmintempfreezedaysout") |
                       ((names(climate)=="Apresepxmintempfreezedaysout") |
              ((names(climate)=="Bresepxmintempfreezedaysout")))))))))))))
if (length(climatehold)!=0) {
  climate <- climate[,-climatehold] # removes the offending columns
}
climate <- data.frame(climate,
                      Aresepxpptout,
                      Apresepxpptout,
                      Bresepxpptout,
                      Aresepxprismjjavpdout,
                      Apresepxprismjjavpdout,
                      Bresepxprismjjavpdout,
                      Aresepxmagnitudeout,
                      Apresepxmagnitudeout,
                      Bresepxmagnitudeout,
                      Aresepxmintempfreezedaysout,
                      Apresepxmintempfreezedaysout,
                      Bresepxmintempfreezedaysout)
#
# Selecting just the resep, resepx and outlier columns
jhold <- which((names(climate)=="Aresep") |
                ((names(climate)=="Apresep") |
                ((names(climate)=="Bresep") |
                ((names(climate)=="Aresepx") |
                ((names(climate)=="Apresepx") |
                ((names(climate)=="Bresepx") |
                ((names(climate)=="Aresepxpptout") |
                ((names(climate)=="Apresepxpptout") |
                ((names(climate)=="Bresepxpptout") |
                ((names(climate)=="Aresepxprismjjavpdout") |
                ((names(climate)=="Apresepxprismjjavpdout") |
                ((names(climate)=="Bresepxprismjjavpdout") |
                ((names(climate)=="Aresepxmagnitudeout") |
                ((names(climate)=="Apresepxmagnitudeout") |
                ((names(climate)=="Bresepxmagnitudeout") |
                ((names(climate)=="Aresepxmintempfreezedaysout") |
                ((names(climate)=="Apresepxmintempfreezedaysout") |
          (names(climate)=="Bresepxmintempfreezedaysout"))))))))))))))))))
climatex <- climate[,c(1:20,jhold)]
#
# Write out the final data for Ruth Kerry to analyze with spatial regression 
# on the aggregated data (climatex) as separate files:
setwd(oname)
#
# A horizons
ihold <- which(!is.na(climatex$Aresepx))
jhold <- which((names(climatex)=="Apresep") |
                 ((names(climatex)=="Bresep") |
                 ((names(climatex)=="Apresepx") |
                 ((names(climatex)=="Bresepx") |
                 ((names(climatex)=="Apresepxpptout") |
                 ((names(climatex)=="Bresepxpptout") |
                 ((names(climatex)=="Apresepxprismjjavpdout") |
                 ((names(climatex)=="Bresepxprismjjavpdout") |
                 ((names(climatex)=="Apresepxmagnitudeout") |
                 ((names(climatex)=="Bresepxmagnitudeout") |
                 ((names(climatex)=="Apresepxmintempfreezedaysout") |
                 (names(climatex)=="Bresepxmintempfreezedaysout"))))))))))))
Aclimatex <- climatex[ihold,-jhold]
names(Aclimatex)[which(names(Aclimatex)=="Aresepx")] <- "resepx"
names(Aclimatex)[which(names(Aclimatex)=="Aresep")] <- "resep"
names(Aclimatex)[
  which(names(Aclimatex)=="Aresepxpptout")] <- 
  "resepxpptout"
names(Aclimatex)[
  which(names(Aclimatex)=="Aresepxprismjjavpdout")] <- 
  "resepxprismjjavpdout"
names(Aclimatex)[
  which(names(Aclimatex)=="Aresepxmagnitudeout")] <- 
  "resepxmagnitudeout"
names(Aclimatex)[
  which(names(Aclimatex)=="Aresepxmintempfreezedaysout")] <- 
  "resepxmintempfreezedaysout"
write.table(Aclimatex,file=paste(gsub("-", "", as.character(Sys.Date())),
            '_EPProject_AHorizons_',
            'AggregatedResultsofSpatiallyModeledRESEPx.csv',sep=""),
            sep=",",row.names=F)
#
# Ap horizons
ihold <- which(!is.na(climatex$Apresepx))
jhold <- which((names(climatex)=="Aresep") |
                ((names(climatex)=="Bresep") |
                ((names(climatex)=="Aresepx") |
                ((names(climatex)=="Bresepx") |
                ((names(climatex)=="Aresepxpptout") |
                ((names(climatex)=="Bresepxpptout") |
                ((names(climatex)=="Aresepxprismjjavpdout") |
                ((names(climatex)=="Bresepxprismjjavpdout") |
                ((names(climatex)=="Aresepxmagnitudeout") |
                ((names(climatex)=="Bresepxmagnitudeout") |
                ((names(climatex)=="Aresepxmintempfreezedaysout") |
                (names(climatex)=="Bresepxmintempfreezedaysout"))))))))))))
Apclimatex <- climatex[ihold,-jhold]
names(Apclimatex)[which(names(Apclimatex)=="Apresepx")] <- "resepx"
names(Apclimatex)[which(names(Apclimatex)=="Apresep")] <- "resep"
names(Apclimatex)[
  which(names(Apclimatex)=="Apresepxpptout")] <- 
  "resepxpptout"
names(Apclimatex)[
  which(names(Apclimatex)=="Apresepxprismjjavpdout")] <- 
  "resepxprismjjavpdout"
names(Apclimatex)[
  which(names(Apclimatex)=="Apresepxmagnitudeout")] <- 
  "resepxmagnitudeout"
names(Apclimatex)[
  which(names(Apclimatex)=="Apresepxmintempfreezedaysout")] <- 
  "resepxmintempfreezedaysout"
write.table(Apclimatex,file=paste(gsub("-", "", as.character(Sys.Date())),
            '_EPProject_ApHorizons_',
            'AggregatedResultsofSpatiallyModeledRESEPx.csv',sep=""),
            sep=",",row.names=F)
#
# B horizons
ihold <- which(!is.na(climatex$Bresepx))
jhold <- which((names(climatex)=="Aresep") |
                ((names(climatex)=="Apresep") |
                ((names(climatex)=="Aresepx") |
                ((names(climatex)=="Apresepx") |
                ((names(climatex)=="Aresepxpptout") |
                ((names(climatex)=="Apresepxpptout") |
                ((names(climatex)=="Aresepxprismjjavpdout") |
                ((names(climatex)=="Apresepxprismjjavpdout") |
                ((names(climatex)=="Aresepxmagnitudeout") |
                ((names(climatex)=="Apresepxmagnitudeout") |
                ((names(climatex)=="Aresepxmintempfreezedaysout") |
                (names(climatex)=="Apresepxmintempfreezedaysout"))))))))))))
Bclimatex <- climatex[ihold,-jhold]
names(Bclimatex)[which(names(Bclimatex)=="Bresepx")] <- "resepx"
names(Bclimatex)[which(names(Bclimatex)=="Bresep")] <- "resep"
names(Bclimatex)[
  which(names(Bclimatex)=="Bresepxpptout")] <- 
  "resepxpptout"
names(Bclimatex)[
  which(names(Bclimatex)=="Bresepxprismjjavpdout")] <- 
  "resepxprismjjavpdout"
names(Bclimatex)[
  which(names(Bclimatex)=="Bresepxmagnitudeout")] <- 
  "resepxmagnitudeout"
names(Bclimatex)[
  which(names(Bclimatex)=="Bresepxmintempfreezedaysout")] <- 
  "resepxmintempfreezedaysout"
write.table(Bclimatex,file=paste(gsub("-", "", as.character(Sys.Date())),
            '_EPProject_BHorizons_',
            'AggregatedResultsofSpatiallyModeledRESEPx.csv',sep=""),
            sep=",",row.names=F)
setwd(fname)
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Read in the new data from Ruth's analysis (3/30/18), select the closest
# sample to each station, sample randomly so that the points are at least a
# certain distance apart (creating at least 5 sets of data), run the 
# regression in Figure 2, and send the 5 sets of data to Ruth along with 
# the results of the analysis (4/28/18)
# ---------------------------------------------------------------------------
# Read in the new data
# --------------------
setwd(dname)
Aclimatex2 <- read.csv("20180509_EPProject_AHorizons_AggregatedResultsofSpatiallyModeledRESEPx.conus.csv") # 
Apclimatex2 <- read.csv("20180509_EPProject_ApHorizons_AggregatedResultsofSpatiallyModeledRESEPx.conus.csv") #
Bclimatex2 <- read.csv("20180509_EPProject_BHorizons_AggregatedResultsofSpatiallyModeledRESEPx.conus.csv") # 
Ahorizons <- read.csv("AHorizons_checked.csv") # 
Aphorizons <- read.csv("ApHorizons_checked.csv") # 
Bhorizons <- read.csv("BHorizons_checked.csv") # 
setwd(wname)
#
# Fixes the comma problem in the station id column
Aclimatex2$stnid <- as.numeric(gsub(",", "", 
                                    as.character(Aclimatex2$stnid)))
Apclimatex2$stnid <- as.numeric(gsub(",", "", 
                                     as.character(Apclimatex2$stnid)))
Bclimatex2$stnid <- as.numeric(gsub(",", "", 
                                    as.character(Bclimatex2$stnid)))
#
#
# Bootstrapping the selection of one random point near each station
# -----------------------------------------------------------------
# Explanation: A subset of points are chosen such that each station with 
# nearby soil samples are paired to exactly one soil sample chosen at 
# random. This set of points is chosen over many interations (bootstrap) 
# to see if the regressions can be run from just those data which because 
# of the work Ruth did earlier with the spatial model, there should not 
# be a significant spatial effect on the regressions. That is, the 
# aggregation appeared to introduce the spatial effect so sampling just 
# one point near the station will remove that aggregation and, thus, 
# spatial effect.
# 
# Examining the data that Ruth sent
map('usa')
map("state", boundary = FALSE, col="gray", add = TRUE)
points(lat~long,Aclimatex2,cex=0.5,pch=21,bg="blue")
points(lat~long,Ahorizons,cex=0.5,pch=22,bg="red")
#
map('usa')
map("state", boundary = FALSE, col="gray", add = TRUE)
points(lat~long,Bclimatex2,cex=0.5,pch=21,bg="blue")
points(lat~long,Bhorizons,cex=0.5,pch=22,bg="red")
#
# The datasets do not seem to be aligned properly here. Something is wrong
# with the Y dimension of the climatex2 and horizons datasets:
plot(Y~X,Apclimatex2,cex=0.5,pch=21,bg="blue",ylim=c(-1500000,3500000))
points(Y~X,Aphorizons,cex=0.5,pch=22,bg="red")
#
# Load the package 'fields' (See note in the preamble above)
require(fields)
#
# Select the nearest station for the A horizons:
#
# Switch the variables to generalize in the code below
myclimatehold <- Aclimatex2
mydhold <- Ahorizons
#
# Create a distance matrix to find the minimum distance
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
da <- as.matrix(data.frame(lon=mydhold$long,lat=mydhold$lat))
distmat <- rdist.earth(cl,da) # rows are the row index for cl; cols are the
#   row index for da
#
#
# Run a loop to record the nearest soil sample to each climate station
stmat <- matrix(NA,nrow=nrow(distmat),ncol=ncol(distmat))
for (j in 1:ncol(distmat)) {
  stmat[which.min(distmat[,j]),j] <- 1
}
#
# Run a loop N times to bootstrap the selection of the nearest points
N <- 1200
myclimatex1random <- matrix(NA,ncol=N,nrow=nrow(myclimatehold)) # here each
# column represents one set of RESEPx values out of all possible sets of 
# RESEPx values closest to each station for all the stations
kb <- txtProgressBar(min = 0, max = N, style = 3) # k progress bar
for (k in 1:N) {
  # Run a loop to choose a residual effective porosity x for the 
  #   A horizons per station
  steff <- NULL
  for (i in 1:nrow(stmat)) {
    hold <- which(stmat[i,]==1)
    if (length(hold) > 1) hold <- sample(hold,1) # randomly selects the 
    # horizon if more than one are present at the pit closest to the station
    if (length(hold)>0) {
      steff <- c(steff,mydhold$RESEPx[hold])
    } else {
      steff <- c(steff,NA)
    }
  }
  myclimatex1random[,k] <- steff
  Sys.sleep(0.1)
  # update GUI console
  setTxtProgressBar(kb, k) 
}
close(kb)
#
Aclimatex1random <- data.frame(myclimatehold,myclimatex1random)
#
#
# Select the nearest station for the Ap horizons:
#
# Switch the variables to generalize in the code below
myclimatehold <- Apclimatex2
mydhold <- Aphorizons
#
# Create a distance matrix to find the minimum distance
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
da <- as.matrix(data.frame(lon=mydhold$long,lat=mydhold$lat))
distmat <- rdist.earth(cl,da) # rows are the row index for cl; cols are the
#   row index for da
#
# Run a loop to record the nearest soil sample to each climate station
stmat <- matrix(NA,nrow=nrow(distmat),ncol=ncol(distmat))
for (j in 1:ncol(distmat)) {
  stmat[which.min(distmat[,j]),j] <- 1
}
#
# Run a loop N times to bootstrap the selection of the nearest points
N <- 1200
myclimatex1random <- matrix(NA,ncol=N,nrow=nrow(myclimatehold)) # here each
# column represents one set of RESEPx values out of all possible sets of 
# RESEPx values closest to each station for all the stations
kb <- txtProgressBar(min = 0, max = N, style = 3) # k progress bar
for (k in 1:N) {
  # Run a loop to choose a residual effective porosity x for the 
  #   Ap horizons per station
  steff <- NULL
  for (i in 1:nrow(stmat)) {
    hold <- which(stmat[i,]==1)
    if (length(hold) > 1) hold <- sample(hold,1) # randomly selects the 
    # horizon if more than one are present at the pit closest to the station
    if (length(hold)>0) {
      steff <- c(steff,mydhold$RESEPx[hold])
    } else {
      steff <- c(steff,NA)
    }
  }
  myclimatex1random[,k] <- steff
  Sys.sleep(0.1)
  # update GUI console
  setTxtProgressBar(kb, k) 
}
close(kb)
#
Apclimatex1random <- data.frame(myclimatehold,myclimatex1random)
#
#
# Select the nearest station for the B horizons:
#
# Switch the variables to generalize in the code below
myclimatehold <- Bclimatex2
mydhold <- Bhorizons
#
# Create a distance matrix to find the minimum distance
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
da <- as.matrix(data.frame(lon=mydhold$long,lat=mydhold$lat))
distmat <- rdist.earth(cl,da) # rows are the row index for cl; cols are the
#   row index for da
#
# Run a loop to record the nearest soil sample to each climate station
stmat <- matrix(NA,nrow=nrow(distmat),ncol=ncol(distmat))
for (j in 1:ncol(distmat)) {
  stmat[which.min(distmat[,j]),j] <- 1
}
#
# Run a loop N times to bootstrap the selection of the nearest points
N <- 1200
myclimatex1random <- matrix(NA,ncol=N,nrow=nrow(myclimatehold)) # here each
# column represents one set of RESEPx values out of all possible sets of 
# RESEPx values closest to each station for all the stations
kb <- txtProgressBar(min = 0, max = N, style = 3) # k progress bar
for (k in 1:N) {
  # Run a loop to choose a residual effective porosity x for the 
  #   B horizons per station
  steff <- NULL
  for (i in 1:nrow(stmat)) {
    hold <- which(stmat[i,]==1)
    if (length(hold) > 1) hold <- sample(hold,1) # randomly selects the 
    # horizon if more than one are present at the pit closest to the station
    if (length(hold)>0) {
      steff <- c(steff,mydhold$RESEPx[hold])
    } else {
      steff <- c(steff,NA)
    }
  }
  myclimatex1random[,k] <- steff
  Sys.sleep(0.1)
  # update GUI console
  setTxtProgressBar(kb, k) 
}
close(kb)
#
Bclimatex1random <- data.frame(myclimatehold,myclimatex1random)
#
#
# Function to extract a pvalue from an lm object
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
#
#
# Run the regressions in Figure 2 with these data (these do not remove the 
# outliers)
#
# Dataframes to store the data
Apptresepx1randomcoef <- data.frame(intercept = rep(NA,N), 
                                    slope = NA,
                                    r2 = NA,
                                    pvalue = NA)
Avpdresepx1randomcoef <- Apptresepx1randomcoef
Amagresepx1randomcoef <- Apptresepx1randomcoef
Afreezeresepx1randomcoef <- Apptresepx1randomcoef
Appptresepx1randomcoef <- Apptresepx1randomcoef
Apvpdresepx1randomcoef <- Apptresepx1randomcoef
Apmagresepx1randomcoef <- Apptresepx1randomcoef
Apfreezeresepx1randomcoef <- Apptresepx1randomcoef
Bpptresepx1randomcoef <- Apptresepx1randomcoef
Bvpdresepx1randomcoef <- Apptresepx1randomcoef
Bmagresepx1randomcoef <- Apptresepx1randomcoef
Bfreezeresepx1randomcoef <- Apptresepx1randomcoef
#
# Select the x data
myAppt <- Aclimatex1random[,"ppt"]
myApppt <- Apclimatex1random[,"ppt"]
myBppt <- Bclimatex1random[,"ppt"]
myAvpd <- Aclimatex1random[,"prismjjavpd"]
myApvpd <- Apclimatex1random[,"prismjjavpd"]
myBvpd <- Bclimatex1random[,"prismjjavpd"]
myAmag <- Aclimatex1random[,"magnitude"]
myApmag <- Apclimatex1random[,"magnitude"]
myBmag <- Bclimatex1random[,"magnitude"]
myAfreeze <- Aclimatex1random[,"mintempfreezedays"]
myApfreeze <- Apclimatex1random[,"mintempfreezedays"]
myBfreeze <- Bclimatex1random[,"mintempfreezedays"]
for (i in 1:N) {
  # Select the y data
  myAresepx <- Aclimatex1random[,paste("X",i,sep="")]
  myApresepx <- Apclimatex1random[,paste("X",i,sep="")]
  myBresepx <- Bclimatex1random[,paste("X",i,sep="")]
  #
  # Run the regression:
  # A horizons
  lmhold <- lm(myAresepx~myAppt)
  Apptresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                 lmp(lmhold))
  lmhold <- lm(myAresepx~myAvpd)
  Avpdresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                 lmp(lmhold))
  lmhold <- lm(myAresepx~myAmag)
  Amagresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                 lmp(lmhold))
  lmhold <- lm(myAresepx~myAfreeze)
  Afreezeresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                    lmp(lmhold))
  #
  # Ap horizons
  lmhold <- lm(myApresepx~myApppt)
  Appptresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                  lmp(lmhold))
  lmhold <- lm(myApresepx~myApvpd)
  Apvpdresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                  lmp(lmhold))
  lmhold <- lm(myApresepx~myApmag)
  Apmagresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                  lmp(lmhold))
  lmhold <- lm(myApresepx~myApfreeze)
  Apfreezeresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                     lmp(lmhold))
  #
  # B horizons
  lmhold <- lm(myBresepx~myBppt)
  Bpptresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                 lmp(lmhold))
  lmhold <- lm(myBresepx~myBvpd)
  Bvpdresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                 lmp(lmhold))
  lmhold <- lm(myBresepx~myBmag)
  Bmagresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                 lmp(lmhold))
  lmhold <- lm(myBresepx~myBfreeze)
  Bfreezeresepx1randomcoef[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                                    lmp(lmhold))
}
#
#
# Aggregate the results of the regression coefficients (these estimates 
# do not currently have outliers removed)
slopemeanmindist1random <- data.frame(horizon=c("A","Ap","B"),
                               ppt=NA,
                               vpd=NA,
                               mag=NA,
                               freeze=NA)
interceptmeanmindist1random <- slopemeanmindist1random
slopesdmindist1random <- slopemeanmindist1random
interceptsdmindist1random <- slopemeanmindist1random
slopePvaluemindist1random <- slopemeanmindist1random
interceptPvaluemindist1random <- slopemeanmindist1random
bootAlist <- list(Apptresepx1randomcoef,
                  Avpdresepx1randomcoef,
                  Amagresepx1randomcoef,
                  Afreezeresepx1randomcoef)
bootAplist <- list(Appptresepx1randomcoef,
                   Apvpdresepx1randomcoef,
                   Apmagresepx1randomcoef,
                   Apfreezeresepx1randomcoef)
bootBlist <- list(Bpptresepx1randomcoef,
                  Bvpdresepx1randomcoef,
                  Bmagresepx1randomcoef,
                  Bfreezeresepx1randomcoef)
for (i in 1:length(bootAlist)) {
  slopemeanmindist1random[1,i+1] <- mean(bootAlist[[i]][,"slope"])
  interceptmeanmindist1random[1,i+1] <- mean(bootAlist[[i]][,"intercept"])
  slopesdmindist1random[1,i+1] <- sd(bootAlist[[i]][,"slope"])
  interceptsdmindist1random[1,i+1] <- sd(bootAlist[[i]][,"intercept"])
  slopePvaluemindist1random[1,i+1] <- t.test(bootAlist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindist1random[1,i+1] <- t.test(bootAlist[[i]][,"intercept"],
                                          mu=0)$p.value
  slopemeanmindist1random[2,i+1] <- mean(bootAplist[[i]][,"slope"])
  interceptmeanmindist1random[2,i+1] <- mean(bootAplist[[i]][,"intercept"])
  slopesdmindist1random[2,i+1] <- sd(bootAplist[[i]][,"slope"])
  interceptsdmindist1random[2,i+1] <- sd(bootAplist[[i]][,"intercept"])
  slopePvaluemindist1random[2,i+1] <- t.test(bootAplist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindist1random[2,i+1] <- t.test(bootAplist[[i]][,"intercept"],
                                          mu=0)$p.value
  slopemeanmindist1random[3,i+1] <- mean(bootBlist[[i]][,"slope"])
  interceptmeanmindist1random[3,i+1] <- mean(bootBlist[[i]][,"intercept"])
  slopesdmindist1random[3,i+1] <- sd(bootBlist[[i]][,"slope"])
  interceptsdmindist1random[3,i+1] <- sd(bootBlist[[i]][,"intercept"])
  slopePvaluemindist1random[3,i+1] <- t.test(bootBlist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindist1random[3,i+1] <- t.test(bootBlist[[i]][,"intercept"],
                                          mu=0)$p.value
}
#
#
# Add outlier information for the X1-X5 columns in a separate dataframe
Aresepx1randompptout <- data.frame(X1out=rep(FALSE,nrow(Aclimatex1random)),
                                   X2out=FALSE,X3out=FALSE,
                                   X4out=FALSE,X5out=FALSE)
Apresepx1randompptout <- data.frame(X1out=rep(FALSE,nrow(Apclimatex1random)),
                                    X2out=FALSE,X3out=FALSE,
                                    X4out=FALSE,X5out=FALSE)
Bresepx1randompptout <- data.frame(X1out=rep(FALSE,nrow(Bclimatex1random)),
                                   X2out=FALSE,X3out=FALSE,
                                   X4out=FALSE,X5out=FALSE)
Aresepx1randomvpdout <- data.frame(X1out=rep(FALSE,nrow(Aclimatex1random)),
                                   X2out=FALSE,X3out=FALSE,
                                   X4out=FALSE,X5out=FALSE)
Apresepx1randomvpdout <- data.frame(X1out=rep(FALSE,nrow(Apclimatex1random)),
                                    X2out=FALSE,X3out=FALSE,
                                    X4out=FALSE,X5out=FALSE)
Bresepx1randomvpdout <- data.frame(X1out=rep(FALSE,nrow(Bclimatex1random)),
                                   X2out=FALSE,X3out=FALSE,
                                   X4out=FALSE,X5out=FALSE)
Aresepx1randommagout <- data.frame(X1out=rep(FALSE,nrow(Aclimatex1random)),
                                   X2out=FALSE,X3out=FALSE,
                                   X4out=FALSE,X5out=FALSE)
Apresepx1randommagout <- data.frame(X1out=rep(FALSE,nrow(Apclimatex1random)),
                                    X2out=FALSE,X3out=FALSE,
                                    X4out=FALSE,X5out=FALSE)
Bresepx1randommagout <- data.frame(X1out=rep(FALSE,nrow(Bclimatex1random)),
                                   X2out=FALSE,X3out=FALSE,
                                   X4out=FALSE,X5out=FALSE)
Aresepx1randomfreezeout <- data.frame(X1out=rep(FALSE,nrow(Aclimatex1random)),
                                      X2out=FALSE,X3out=FALSE,
                                      X4out=FALSE,X5out=FALSE)
Apresepx1randomfreezeout <- data.frame(X1out=rep(FALSE,nrow(Apclimatex1random)),
                                       X2out=FALSE,X3out=FALSE,
                                       X4out=FALSE,X5out=FALSE)
Bresepx1randomfreezeout <- data.frame(X1out=rep(FALSE,nrow(Bclimatex1random)),
                                      X2out=FALSE,X3out=FALSE,
                                      X4out=FALSE,X5out=FALSE)
#
#
# Function to identify bivariate outliers using the adjusted Mahalanobois 
# distance procedure
outlierID <- function(df) {
  # Creates a logical vector to test to see if the point is a bivariate
  # outlier; vectore corresponds to the order of the rows in the input df
  # df = data.frame; only the first two columns will be used to detect the
  # outliers with; those first two columns need to be numeric
  #
  # Please note: there is a bizzare problem with mvOutlier where if the 
  # row.names of the dataframe are NULL from 1 to nrow(df), the algorithm
  # will designated the number of outliers correctly but will remove the
  # wrong ones in the resulting dfout$newData dataframe. Solution is to 
  row.names(df) <- 1:nrow(df)
  dfout <- mvOutlier(df,qqplot=F,method="adj.quan")
  dfo <- dfout$outlier # data.frame indicating the outliers
  dfnd <- dfout$newData # df without the outliers
  ifcheck <- length(which(dfo$Outlier==TRUE)) > 0
  if (ifcheck) {
    out <- rep(FALSE,nrow(df))
    dfindex <- as.numeric(row.names(dfo[which(dfo$Outlier==TRUE),]))
    out[dfindex] <- TRUE
  } else {
    out <- rep(FALSE,nrow(df))
  }
  out # TRUE indicates outlier; FALSE indicates not an outlier
}
#
#
# MAP vs resepx for A horizons
x <- Aclimatex1random$ppt
resepxX1 <- Aclimatex1random$X1
resepxX2 <- Aclimatex1random$X2
resepxX3 <- Aclimatex1random$X3
resepxX4 <- Aclimatex1random$X4
resepxX5 <- Aclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Aresepx1randompptout$X1out <- outlierID(dholdX1)
Aresepx1randompptout$X2out <- outlierID(dholdX2)
Aresepx1randompptout$X3out <- outlierID(dholdX3)
Aresepx1randompptout$X4out <- outlierID(dholdX4)
Aresepx1randompptout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~ppt,Aclimatex1random[which(Aresepx1randompptout$X1out==FALSE),])
plot(X2~ppt,Aclimatex1random[which(Aresepx1randompptout$X2out==FALSE),])
plot(X3~ppt,Aclimatex1random[which(Aresepx1randompptout$X3out==FALSE),])
plot(X4~ppt,Aclimatex1random[which(Aresepx1randompptout$X4out==FALSE),])
plot(X5~ppt,Aclimatex1random[which(Aresepx1randompptout$X5out==FALSE),])
#
# MAP vs resepx for Ap horizons
x <- Apclimatex1random$ppt
resepxX1 <- Apclimatex1random$X1
resepxX2 <- Apclimatex1random$X2
resepxX3 <- Apclimatex1random$X3
resepxX4 <- Apclimatex1random$X4
resepxX5 <- Apclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Apresepx1randompptout$X1out <- outlierID(dholdX1)
Apresepx1randompptout$X2out <- outlierID(dholdX2)
Apresepx1randompptout$X3out <- outlierID(dholdX3)
Apresepx1randompptout$X4out <- outlierID(dholdX4)
Apresepx1randompptout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~ppt,Apclimatex1random[which(Apresepx1randompptout$X1out==FALSE),])
plot(X2~ppt,Apclimatex1random[which(Apresepx1randompptout$X2out==FALSE),])
plot(X3~ppt,Apclimatex1random[which(Apresepx1randompptout$X3out==FALSE),])
plot(X4~ppt,Apclimatex1random[which(Apresepx1randompptout$X4out==FALSE),])
plot(X5~ppt,Apclimatex1random[which(Apresepx1randompptout$X5out==FALSE),])
#
# MAP vs resepx for B horizons
x <- Bclimatex1random$ppt
resepxX1 <- Bclimatex1random$X1
resepxX2 <- Bclimatex1random$X2
resepxX3 <- Bclimatex1random$X3
resepxX4 <- Bclimatex1random$X4
resepxX5 <- Bclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Bresepx1randompptout$X1out <- outlierID(dholdX1)
Bresepx1randompptout$X2out <- outlierID(dholdX2)
Bresepx1randompptout$X3out <- outlierID(dholdX3)
Bresepx1randompptout$X4out <- outlierID(dholdX4)
Bresepx1randompptout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~ppt,Bclimatex1random[which(Bresepx1randompptout$X1out==FALSE),])
plot(X2~ppt,Bclimatex1random[which(Bresepx1randompptout$X2out==FALSE),])
plot(X3~ppt,Bclimatex1random[which(Bresepx1randompptout$X3out==FALSE),])
plot(X4~ppt,Bclimatex1random[which(Bresepx1randompptout$X4out==FALSE),])
plot(X5~ppt,Bclimatex1random[which(Bresepx1randompptout$X5out==FALSE),])
#
# VPD vs resepx for A horizons
x <- Aclimatex1random$prismjjavpd
resepxX1 <- Aclimatex1random$X1
resepxX2 <- Aclimatex1random$X2
resepxX3 <- Aclimatex1random$X3
resepxX4 <- Aclimatex1random$X4
resepxX5 <- Aclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Aresepx1randomvpdout$X1out <- outlierID(dholdX1)
Aresepx1randomvpdout$X2out <- outlierID(dholdX2)
Aresepx1randomvpdout$X3out <- outlierID(dholdX3)
Aresepx1randomvpdout$X4out <- outlierID(dholdX4)
Aresepx1randomvpdout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~prismjjavpd,Aclimatex1random[which(Aresepx1randomvpdout$X1out==FALSE),])
plot(X2~prismjjavpd,Aclimatex1random[which(Aresepx1randomvpdout$X2out==FALSE),])
plot(X3~prismjjavpd,Aclimatex1random[which(Aresepx1randomvpdout$X3out==FALSE),])
plot(X4~prismjjavpd,Aclimatex1random[which(Aresepx1randomvpdout$X4out==FALSE),])
plot(X5~prismjjavpd,Aclimatex1random[which(Aresepx1randomvpdout$X5out==FALSE),])
#
# VPD vs resepx for Ap horizons
x <- Apclimatex1random$prismjjavpd
resepxX1 <- Apclimatex1random$X1
resepxX2 <- Apclimatex1random$X2
resepxX3 <- Apclimatex1random$X3
resepxX4 <- Apclimatex1random$X4
resepxX5 <- Apclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Apresepx1randomvpdout$X1out <- outlierID(dholdX1)
Apresepx1randomvpdout$X2out <- outlierID(dholdX2)
Apresepx1randomvpdout$X3out <- outlierID(dholdX3)
Apresepx1randomvpdout$X4out <- outlierID(dholdX4)
Apresepx1randomvpdout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~prismjjavpd,Apclimatex1random[which(Apresepx1randomvpdout$X1out==FALSE),])
plot(X2~prismjjavpd,Apclimatex1random[which(Apresepx1randomvpdout$X2out==FALSE),])
plot(X3~prismjjavpd,Apclimatex1random[which(Apresepx1randomvpdout$X3out==FALSE),])
plot(X4~prismjjavpd,Apclimatex1random[which(Apresepx1randomvpdout$X4out==FALSE),])
plot(X5~prismjjavpd,Apclimatex1random[which(Apresepx1randomvpdout$X5out==FALSE),])
#
# VPD vs resepx for B horizons
x <- Bclimatex1random$prismjjavpd
resepxX1 <- Bclimatex1random$X1
resepxX2 <- Bclimatex1random$X2
resepxX3 <- Bclimatex1random$X3
resepxX4 <- Bclimatex1random$X4
resepxX5 <- Bclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Bresepx1randomvpdout$X1out <- outlierID(dholdX1)
Bresepx1randomvpdout$X2out <- outlierID(dholdX2)
Bresepx1randomvpdout$X3out <- outlierID(dholdX3)
Bresepx1randomvpdout$X4out <- outlierID(dholdX4)
Bresepx1randomvpdout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~prismjjavpd,Bclimatex1random[which(Bresepx1randomvpdout$X1out==FALSE),])
plot(X2~prismjjavpd,Bclimatex1random[which(Bresepx1randomvpdout$X2out==FALSE),])
plot(X3~prismjjavpd,Bclimatex1random[which(Bresepx1randomvpdout$X3out==FALSE),])
plot(X4~prismjjavpd,Bclimatex1random[which(Bresepx1randomvpdout$X4out==FALSE),])
plot(X5~prismjjavpd,Bclimatex1random[which(Bresepx1randomvpdout$X5out==FALSE),])
#
# Precip magnitude vs resepx for A horizons
x <- Aclimatex1random$magnitude
resepxX1 <- Aclimatex1random$X1
resepxX2 <- Aclimatex1random$X2
resepxX3 <- Aclimatex1random$X3
resepxX4 <- Aclimatex1random$X4
resepxX5 <- Aclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Aresepx1randommagout$X1out <- outlierID(dholdX1)
Aresepx1randommagout$X2out <- outlierID(dholdX2)
Aresepx1randommagout$X3out <- outlierID(dholdX3)
Aresepx1randommagout$X4out <- outlierID(dholdX4)
Aresepx1randommagout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~magnitude,Aclimatex1random[which(Aresepx1randommagout$X1out==FALSE),])
plot(X2~magnitude,Aclimatex1random[which(Aresepx1randommagout$X2out==FALSE),])
plot(X3~magnitude,Aclimatex1random[which(Aresepx1randommagout$X3out==FALSE),])
plot(X4~magnitude,Aclimatex1random[which(Aresepx1randommagout$X4out==FALSE),])
plot(X5~magnitude,Aclimatex1random[which(Aresepx1randommagout$X5out==FALSE),])
#
# Precip magnitude vs resepx for Ap horizons
x <- Apclimatex1random$magnitude
resepxX1 <- Apclimatex1random$X1
resepxX2 <- Apclimatex1random$X2
resepxX3 <- Apclimatex1random$X3
resepxX4 <- Apclimatex1random$X4
resepxX5 <- Apclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Apresepx1randommagout$X1out <- outlierID(dholdX1)
Apresepx1randommagout$X2out <- outlierID(dholdX2)
Apresepx1randommagout$X3out <- outlierID(dholdX3)
Apresepx1randommagout$X4out <- outlierID(dholdX4)
Apresepx1randommagout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~magnitude,Apclimatex1random[which(Apresepx1randommagout$X1out==FALSE),])
plot(X2~magnitude,Apclimatex1random[which(Apresepx1randommagout$X2out==FALSE),])
plot(X3~magnitude,Apclimatex1random[which(Apresepx1randommagout$X3out==FALSE),])
plot(X4~magnitude,Apclimatex1random[which(Apresepx1randommagout$X4out==FALSE),])
plot(X5~magnitude,Apclimatex1random[which(Apresepx1randommagout$X5out==FALSE),])
#
# Precip magnitude vs resepx for B horizons
x <- Bclimatex1random$magnitude
resepxX1 <- Bclimatex1random$X1
resepxX2 <- Bclimatex1random$X2
resepxX3 <- Bclimatex1random$X3
resepxX4 <- Bclimatex1random$X4
resepxX5 <- Bclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Bresepx1randommagout$X1out <- outlierID(dholdX1)
Bresepx1randommagout$X2out <- outlierID(dholdX2)
Bresepx1randommagout$X3out <- outlierID(dholdX3)
Bresepx1randommagout$X4out <- outlierID(dholdX4)
Bresepx1randommagout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~magnitude,Bclimatex1random[which(Bresepx1randommagout$X1out==FALSE),])
plot(X2~magnitude,Bclimatex1random[which(Bresepx1randommagout$X2out==FALSE),])
plot(X3~magnitude,Bclimatex1random[which(Bresepx1randommagout$X3out==FALSE),])
plot(X4~magnitude,Bclimatex1random[which(Bresepx1randommagout$X4out==FALSE),])
plot(X5~magnitude,Bclimatex1random[which(Bresepx1randommagout$X5out==FALSE),])
#
# Freeze frequency vs resepx for A horizons
x <- Aclimatex1random$mintempfreezedays
resepxX1 <- Aclimatex1random$X1
resepxX2 <- Aclimatex1random$X2
resepxX3 <- Aclimatex1random$X3
resepxX4 <- Aclimatex1random$X4
resepxX5 <- Aclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Aresepx1randomfreezeout$X1out <- outlierID(dholdX1)
Aresepx1randomfreezeout$X2out <- outlierID(dholdX2)
Aresepx1randomfreezeout$X3out <- outlierID(dholdX3)
Aresepx1randomfreezeout$X4out <- outlierID(dholdX4)
Aresepx1randomfreezeout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~mintempfreezedays,Aclimatex1random[which(
  Aresepx1randomfreezeout$X1out==FALSE),])
plot(X2~mintempfreezedays,Aclimatex1random[which(
  Aresepx1randomfreezeout$X2out==FALSE),])
plot(X3~mintempfreezedays,Aclimatex1random[which(
  Aresepx1randomfreezeout$X3out==FALSE),])
plot(X4~mintempfreezedays,Aclimatex1random[which(
  Aresepx1randomfreezeout$X4out==FALSE),])
plot(X5~mintempfreezedays,Aclimatex1random[which(
  Aresepx1randomfreezeout$X5out==FALSE),])
#
# Freeze frequency vs resepx for Ap horizons
x <- Apclimatex1random$mintempfreezedays
resepxX1 <- Apclimatex1random$X1
resepxX2 <- Apclimatex1random$X2
resepxX3 <- Apclimatex1random$X3
resepxX4 <- Apclimatex1random$X4
resepxX5 <- Apclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Apresepx1randomfreezeout$X1out <- outlierID(dholdX1)
Apresepx1randomfreezeout$X2out <- outlierID(dholdX2)
Apresepx1randomfreezeout$X3out <- outlierID(dholdX3)
Apresepx1randomfreezeout$X4out <- outlierID(dholdX4)
Apresepx1randomfreezeout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~mintempfreezedays,Apclimatex1random[which(
  Apresepx1randomfreezeout$X1out==FALSE),])
plot(X2~mintempfreezedays,Apclimatex1random[which(
  Apresepx1randomfreezeout$X2out==FALSE),])
plot(X3~mintempfreezedays,Apclimatex1random[which(
  Apresepx1randomfreezeout$X3out==FALSE),])
plot(X4~mintempfreezedays,Apclimatex1random[which(
  Apresepx1randomfreezeout$X4out==FALSE),])
plot(X5~mintempfreezedays,Apclimatex1random[which(
  Apresepx1randomfreezeout$X5out==FALSE),])
#
# Freeze frequency vs resepx for B horizons
x <- Bclimatex1random$mintempfreezedays
resepxX1 <- Bclimatex1random$X1
resepxX2 <- Bclimatex1random$X2
resepxX3 <- Bclimatex1random$X3
resepxX4 <- Bclimatex1random$X4
resepxX5 <- Bclimatex1random$X5
dholdX1 <- data.frame(x=x,resepx=resepxX1)
dholdX2 <- data.frame(x=x,resepx=resepxX2)
dholdX3 <- data.frame(x=x,resepx=resepxX3)
dholdX4 <- data.frame(x=x,resepx=resepxX4)
dholdX5 <- data.frame(x=x,resepx=resepxX5)
Bresepx1randomfreezeout$X1out <- outlierID(dholdX1)
Bresepx1randomfreezeout$X2out <- outlierID(dholdX2)
Bresepx1randomfreezeout$X3out <- outlierID(dholdX3)
Bresepx1randomfreezeout$X4out <- outlierID(dholdX4)
Bresepx1randomfreezeout$X5out <- outlierID(dholdX5)
#
# Check
par(mfrow=c(2,4))
plot(X1~mintempfreezedays,Bclimatex1random[which(
  Bresepx1randomfreezeout$X1out==FALSE),])
plot(X2~mintempfreezedays,Bclimatex1random[which(
  Bresepx1randomfreezeout$X2out==FALSE),])
plot(X3~mintempfreezedays,Bclimatex1random[which(
  Bresepx1randomfreezeout$X3out==FALSE),])
plot(X4~mintempfreezedays,Bclimatex1random[which(
  Bresepx1randomfreezeout$X4out==FALSE),])
plot(X5~mintempfreezedays,Bclimatex1random[which(
  Bresepx1randomfreezeout$X5out==FALSE),])
#
#
#
# Rename the columns
names(Aresepx1randompptout) <- paste("ppt",names(Apresepx1randompptout), sep="")
names(Aresepx1randomvpdout) <- paste("vpd",names(Apresepx1randomvpdout), sep="")
names(Aresepx1randommagout) <- paste("mag",names(Apresepx1randommagout), sep="")
names(Aresepx1randomfreezeout) <- paste("freeze",names(Apresepx1randomfreezeout), 
                                        sep="")
names(Apresepx1randompptout) <- paste("ppt",names(Apresepx1randompptout), sep="")
names(Apresepx1randomvpdout) <- paste("vpd",names(Apresepx1randomvpdout), sep="")
names(Apresepx1randommagout) <- paste("mag",names(Apresepx1randommagout), sep="")
names(Apresepx1randomfreezeout) <- paste("freeze",names(Apresepx1randomfreezeout), 
                                        sep="")
names(Bresepx1randompptout) <- paste("ppt",names(Bresepx1randompptout), sep="")
names(Bresepx1randomvpdout) <- paste("vpd",names(Bresepx1randomvpdout), sep="")
names(Bresepx1randommagout) <- paste("mag",names(Bresepx1randommagout), sep="")
names(Bresepx1randomfreezeout) <- paste("freeze",names(Bresepx1randomfreezeout), 
                                        sep="")
#
#
# Output the new files that represent a single random RESEPx near each 
# station
setwd(oname)
# Summary data
write.table(slopemeanmindist1random,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_slopemean1random.csv',sep=""),
            sep=",",row.names=F)
write.table(interceptmeanmindist1random,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_interceptmean1random.csv',sep=""),
            sep=",",row.names=F)
write.table(slopesdmindist1random,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_slopesd1random.csv',sep=""),
            sep=",",row.names=F)
write.table(interceptsdmindist1random,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_interceptsd1random.csv',sep=""),
            sep=",",row.names=F)
write.table(slopePvaluemindist1random,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_slopePvalue1random.csv',sep=""),
            sep=",",row.names=F)

write.table(interceptPvaluemindist1random,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_interceptPvalue1random.csv',sep=""),
            sep=",",row.names=F)
#
# Full bootstrapped data
write.table(Aclimatex1random,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatex1random.csv',sep=""),
            sep=",",row.names=F)
write.table(Apclimatex1random,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatex1random.csv',sep=""),
            sep=",",row.names=F)
write.table(Bclimatex1random,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatex1random.csv',sep=""),
            sep=",",row.names=F)
#
# First 5 bootstrapped data for Ruth
write.table(data.frame(Aclimatex1random[,1:37],
                       Aresepx1randompptout,
                       Aresepx1randomvpdout,
                       Aresepx1randommagout,
                       Aresepx1randomfreezeout),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatex1random5.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Apclimatex1random[,1:37],
                       Apresepx1randompptout,
                       Apresepx1randomvpdout,
                       Apresepx1randommagout,
                       Apresepx1randomfreezeout),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatex1random5.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Bclimatex1random[,1:37],
                       Bresepx1randompptout,
                       Bresepx1randomvpdout,
                       Bresepx1randommagout,
                       Bresepx1randomfreezeout),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatex1random5.csv',sep=""),
            sep=",",row.names=F)
#
setwd(wname)
# 
#
#
# Selecting the closest point to the station
# ------------------------------------------
# Explanation: The subset of points closest to the station will be used
# to see if the regressions can be run from just those data which because
# of the work Ruth did earlier with the spatial model, there should not be
# a significant spatial effect on the regressions. That is, the aggregation
# appeared to introduce the spatial effect so sampling just one point
# near the station will remove that aggregation effect and thus spatial
# effect.
#
# Select an A horizon point for each station:
#
# Switch the variables to generalize in the code below
myclimatehold <- Aclimatex2
mydhold <- Ahorizons
#
# Create a distance matrix to find the minimum distance
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
da <- as.matrix(data.frame(lon=mydhold$long,lat=mydhold$lat))
distmat <- rdist.earth(cl,da) # rows are the row index for cl; cols are the
# row index for da
#
# Run a loop to record the nearest soil sample to each climate station
stmat <- matrix(NA,nrow=nrow(distmat),ncol=ncol(distmat))
for (j in 1:ncol(distmat)) {
  stmat[which.min(distmat[,j]),j] <- 1
}
#
# Run a loop to choose one residual effective porosity x for the 
#   A horizons per station
steff <- NULL
for (i in 1:nrow(stmat)) {
  # Selects the closest sample to the station
  hold <- which(stmat[i,]==1)
  hold <- hold[which(distmat[i,which(stmat[i,]==1)] == 
                       min(distmat[i,which(stmat[i,]==1)]))]
  if (length(hold) > 1) hold <- sample(hold,1) # randomly selects the 
  # horizon if more than one are present at the pit closest to the station
  if (length(hold)>0) {
    steff <- c(steff,mydhold$RESEPx[hold])
  } else {
    steff <- c(steff,NA)
  }
}
myAclimatex1closest <- data.frame(myclimatehold,
                                  closeststatAresepx=steff)
#
#
# Select an Ap horizon point for each station:
#
# Switch the variables to generalize in the code below
myclimatehold <- Apclimatex2
mydhold <- Aphorizons
#
# Create a distance matrix to find the minimum distance
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
da <- as.matrix(data.frame(lon=mydhold$long,lat=mydhold$lat))
distmat <- rdist.earth(cl,da) # rows are the row index for cl; cols are the
# row index for da
#
# Run a loop to record the nearest soil sample to each climate station
stmat <- matrix(NA,nrow=nrow(distmat),ncol=ncol(distmat))
for (j in 1:ncol(distmat)) {
  stmat[which.min(distmat[,j]),j] <- 1
}
#
# Run a loop to choose one residual effective porosity x for the 
#   Ap horizons per station
steff <- NULL
for (i in 1:nrow(stmat)) {
  # Selects the closest sample to the station
  hold <- which(stmat[i,]==1)
  hold <- hold[which(distmat[i,which(stmat[i,]==1)] == 
                       min(distmat[i,which(stmat[i,]==1)]))]
  if (length(hold) > 1) hold <- sample(hold,1) # randomly selects the 
  # horizon if more than one are present at the pit closest to the station
  if (length(hold)>0) {
    steff <- c(steff,mydhold$RESEPx[hold])
  } else {
    steff <- c(steff,NA)
  }
}
myApclimatex1closest <- data.frame(myclimatehold,
                                   closeststatApresepx=steff)
#
#
# Select an B horizon point for each station:
#
# Switch the variables to generalize in the code below
myclimatehold <- Bclimatex2
mydhold <- Bhorizons
#
# Create a distance matrix to find the minimum distance
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
da <- as.matrix(data.frame(lon=mydhold$long,lat=mydhold$lat))
distmat <- rdist.earth(cl,da) # rows are the row index for cl; cols are the
# row index for da
#
# Run a loop to record the nearest soil sample to each climate station
stmat <- matrix(NA,nrow=nrow(distmat),ncol=ncol(distmat))
for (j in 1:ncol(distmat)) {
  stmat[which.min(distmat[,j]),j] <- 1
}
#
# Run a loop to choose one residual effective porosity x for the 
#   B horizons per station
steff <- NULL
for (i in 1:nrow(stmat)) {
  # Selects the closest sample to the station
  hold <- which(stmat[i,]==1)
  hold <- hold[which(distmat[i,which(stmat[i,]==1)] == 
                       min(distmat[i,which(stmat[i,]==1)]))]
  if (length(hold) > 1) hold <- sample(hold,1) # randomly selects the 
  # horizon if more than one are present at the pit closest to the station
  if (length(hold)>0) {
    steff <- c(steff,mydhold$RESEPx[hold])
  } else {
    steff <- c(steff,NA)
  }
}
myBclimatex1closest <- data.frame(myclimatehold,
                                  closeststatBresepx=steff)
#
#
# Detects and identifies outliers in the data
# Add outlier information for the X1 closest columns in a separate dataframe
Aresepx1closestpptout <- data.frame(X1out=rep(FALSE,nrow(myAclimatex1closest)))
Apresepx1closestpptout <- data.frame(X1out=rep(FALSE,nrow(myApclimatex1closest)))
Bresepx1closestpptout <- data.frame(X1out=rep(FALSE,nrow(myBclimatex1closest)))
Aresepx1closestvpdout <- data.frame(X1out=rep(FALSE,nrow(myAclimatex1closest)))
Apresepx1closestvpdout <- data.frame(X1out=rep(FALSE,nrow(myApclimatex1closest)))
Bresepx1closestvpdout <- data.frame(X1out=rep(FALSE,nrow(myBclimatex1closest)))
Aresepx1closestmagout <- data.frame(X1out=rep(FALSE,nrow(myAclimatex1closest)))
Apresepx1closestmagout <- data.frame(X1out=rep(FALSE,nrow(myApclimatex1closest)))
Bresepx1closestmagout <- data.frame(X1out=rep(FALSE,nrow(myBclimatex1closest)))
Aresepx1closestfreezeout <- data.frame(X1out=rep(FALSE,nrow(myAclimatex1closest)))
Apresepx1closestfreezeout <- data.frame(X1out=rep(FALSE,nrow(myApclimatex1closest)))
Bresepx1closestfreezeout <- data.frame(X1out=rep(FALSE,nrow(myBclimatex1closest)))
#
# Rename the resepx variable:
names(myAclimatex1closest)[which(names(myAclimatex1closest)==
                                   "closeststatAresepx")] <- "X1"
names(myApclimatex1closest)[which(names(myApclimatex1closest)==
                                   "closeststatApresepx")] <- "X1"
names(myBclimatex1closest)[which(names(myBclimatex1closest)==
                                   "closeststatBresepx")] <- "X1"
#
# MAP vs resepx for A horizons
x <- myAclimatex1closest$ppt
resepxX1 <- myAclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Aresepx1closestpptout$X1out <- outlierID(dholdX1)
#
# MAP vs resepx for Ap horizons
x <- myApclimatex1closest$ppt
resepxX1 <- myApclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Apresepx1closestpptout$X1out <- outlierID(dholdX1)
#
# MAP vs resepx for B horizons
x <- myBclimatex1closest$ppt
resepxX1 <- myBclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Bresepx1closestpptout$X1out <- outlierID(dholdX1)
#
#
# VPD vs resepx for A horizons
x <- myAclimatex1closest$prismjjavpd
resepxX1 <- myAclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Aresepx1closestvpdout$X1out <- outlierID(dholdX1)
#
# VPD vs resepx for Ap horizons
x <- myApclimatex1closest$prismjjavpd
resepxX1 <- myApclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Apresepx1closestvpdout$X1out <- outlierID(dholdX1)
#
# VPD vs resepx for B horizons
x <- myBclimatex1closest$prismjjavpd
resepxX1 <- myBclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Bresepx1closestvpdout$X1out <- outlierID(dholdX1)
#
# Precip magnitude vs resepx for A horizons
x <- myAclimatex1closest$magnitude
resepxX1 <- myAclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Aresepx1closestmagout$X1out <- outlierID(dholdX1)
#
# Precip magnitude vs resepx for Ap horizons
x <- myApclimatex1closest$magnitude
resepxX1 <- myApclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Apresepx1closestmagout$X1out <- outlierID(dholdX1)
#
# Precip magnitude vs resepx for B horizons
x <- myBclimatex1closest$magnitude
resepxX1 <- myBclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Bresepx1closestmagout$X1out <- outlierID(dholdX1)
#
# Freeze frequency vs resepx for A horizons
x <- myAclimatex1closest$mintempfreezedays
resepxX1 <- myAclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Aresepx1closestfreezeout$X1out <- outlierID(dholdX1)
#
# Freeze frequency vs resepx for Ap horizons
x <- myApclimatex1closest$mintempfreezedays
resepxX1 <- myApclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Apresepx1closestfreezeout$X1out <- outlierID(dholdX1)
#
# Freeze frequency vs resepx for B horizons
x <- myBclimatex1closest$mintempfreezedays
resepxX1 <- myBclimatex1closest$X1
dholdX1 <- data.frame(x=x,resepx=resepxX1)
Bresepx1closestfreezeout$X1out <- outlierID(dholdX1)
#
# Rename the columns
names(Aresepx1closestpptout) <- paste("ppt",names(Apresepx1closestpptout), sep="")
names(Aresepx1closestvpdout) <- paste("vpd",names(Apresepx1closestvpdout), sep="")
names(Aresepx1closestmagout) <- paste("mag",names(Apresepx1closestmagout), sep="")
names(Aresepx1closestfreezeout) <- paste("freeze",names(Apresepx1closestfreezeout), 
                                        sep="")
names(Apresepx1closestpptout) <- paste("ppt",names(Apresepx1closestpptout), sep="")
names(Apresepx1closestvpdout) <- paste("vpd",names(Apresepx1closestvpdout), sep="")
names(Apresepx1closestmagout) <- paste("mag",names(Apresepx1closestmagout), sep="")
names(Apresepx1closestfreezeout) <- paste("freeze",names(Apresepx1closestfreezeout), 
                                         sep="")
names(Bresepx1closestpptout) <- paste("ppt",names(Bresepx1closestpptout), sep="")
names(Bresepx1closestvpdout) <- paste("vpd",names(Bresepx1closestvpdout), sep="")
names(Bresepx1closestmagout) <- paste("mag",names(Bresepx1closestmagout), sep="")
names(Bresepx1closestfreezeout) <- paste("freeze",names(Bresepx1closestfreezeout), 
                                        sep="")
#
#
# Run the regressions in Figure 2 with these data
# [insert code here]
#
#
# Output the three new files that represent a single RESEPx from the closest 
# station.
setwd(oname)
# Closest soil sample to the station for Ruth
write.table(data.frame(myAclimatex1closest,
                       Aresepx1closestpptout,
                       Aresepx1closestvpdout,
                       Aresepx1closestmagout,
                       Aresepx1closestfreezeout),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatex1closest.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(myApclimatex1closest,
                       Apresepx1closestpptout,
                       Apresepx1closestvpdout,
                       Apresepx1closestmagout,
                       Apresepx1closestfreezeout),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatex1closest.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(myBclimatex1closest,
                       Bresepx1closestpptout,
                       Bresepx1closestvpdout,
                       Bresepx1closestmagout,
                       Bresepx1closestfreezeout),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatex1closest.csv',sep=""),
            sep=",",row.names=F)
#
setwd(wname)
#
#
#
# Randomly sample the data so that samples are separated by some minimum 
# distance
# ----------------------------------------------------------------------
# Explanation: The minimum distances were set by Ruth's calculations in
# the previous analysis. Samples must be separated by at least 600 km
# for A horizons, 220 km for Ap horizons, and 350 km for B horizons. This
# is run using the already aggregated data and sampled many times to 
# estimate the parameters of the regression. A climate station with 
# aggregated data is chosen at random and the order of theothers are chosen 
# random in each iteration so that the minimum distancesfollow the results 
# in Ruth's analysis. Each realization of the subset with the appropriate 
# distance is used to estimate the regression parameters.
#
# What follows is large chunks of code for each horizon with the data
# outputted at the very end of the analysis.
#
# Note: there are several functions defined in the A horizon code below
# so those functions are not redefined in the Ap or B horizon code below
# that.
#
# Note: these have been updated (5/9/18) to 300, 100, and 420 km for
# A, Ap, and B horizons respectively; PLEASE NOTE!!!: Naming of the 
# data.frame follows the previous range values even though the ranges
# have been updated with Ruth's analysis on 5/9/18
#
# Subsetting the aggregated A horizon RESEPx data with a minimum distance
# -----------------------------------------------------------------------
#
# Switch the variables to generalize in the code below
mindistance <- 300 # km; for aggregated A horizons
mindistance <- mindistance*1000 # m
myclimatehold <- Aclimatex2
#
# Need the function distHaversine from the  package
require(geosphere)
#
# Create a distance matrix to find the minimum distance in meters
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
distmat <- matrix(NA, nrow=nrow(cl),ncol=nrow(cl))
for (j in 1:nrow(cl)) { # columns in distmat
  for (i in 1:nrow(cl)) { # rows in distmat
    distmat[i,j] <- distHaversine(cl[j,],cl[i,]) # calculates the distance
    # on the ground in meters using the Haversine formula
  }
}
#
#
# Function for finding the indices in a vector that correspond to multiple 
# conditions:
ms <- function(v,x,zequal="equal") {
  # v is a vector to search to see if the conditions are met
  # x is a vector of the conditions
  # zequal is either "equal", "greaterthan", or "lessthan"
  if (zequal == "equal") {
    out <- NULL
    for (i in 1:length(x)) {
      out <- c(out,which(v == x[i]))
    }
    out
  } else {
    if (zequal == "greaterthan") {
      out <- NULL
      for (i in 1:length(x)) {
        out <- c(out,which(v > x[i]))
      }
      out
    } else { # lessthan
      out <- NULL
      for (i in 1:length(x)) {
        out <- c(out,which(v < x[i]))
      }
      out
    }
  }
}
#
#
# Run a loop to randomly select climate stations that are separated by
# a minimum distance. The data for these stations are the aggregated A 
# horizon data after outliers were removed prior to Ruth's previous
# analyses. The loop randomly selects these climate stations by using
# different starting points and order of stations checked each time.
# The loop is run of N times following the code below:
#
N <- 1200 # recommended by Lander as a good number for a general
# bootstrap
clstatlist <- list(NULL) # will hold the results of the loop below
# (run N times) in each element; results are the distmat indices
# that are at the appropriate spacing
#
for (i in 1:N) {
  # Run a loop to select climate stations with a minimum distance
  rdistmatcolindex <- sample(1:ncol(distmat)) # randomizes the column index
  # for distmat so that it doesn't keep starting in the same place
  myindexhold <- sample(which(distmat[,rdistmatcolindex[1]] > mindistance))
  # forms the initial maximum subset of the data
  for (j in 2:length(myindexhold)) {
    mycheck <- which(distmat[,myindexhold[j]] < 
                       mindistance) # holds the distmat indices of the 
    # stations that are too close to the jth myindexhold station under 
    # consideration
    mycheck <- mycheck[-which(mycheck == myindexhold[j])] # removes the
    # distmat index of the jth myindexhold station under consideration
    myhold <- NULL # resets myhold
    #
    myhold <- ms(myindexhold,mycheck) # holds the indices of the 
    # myindexhold vector that should be removed
    #
    if (length(myhold) == 0) myhold <- NULL # no indices were found; 
    # reset myhold to NULL for the if statement below
    #
    if(!is.null(myhold)) myindexhold <- myindexhold[-myhold]
    if (j >= length(myindexhold)) break
  }
  clstatlist[[i]] <- myindexhold
}
#
#
# Check the distribution of all the indices to see how many of the climate
# stations are included in the analysis (bootstrapped so should be most of them):
map('usa')
map("state", boundary = FALSE, col="gray", add = TRUE)
mycolors <- jet.colors(length(clstatlist))
for (i in 1:length(clstatlist)){
  points(lat~long,Aclimatex2[clstatlist[[i]],],cex=0.5,pch=21,bg=mycolors[i])
}
#
#
# Function to extract a pvalue from an lm object
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
#
#
# Run the regressions for each of the appropriately-spaced subsets in
# clstatlist:
bootpptcoefhold <- data.frame(intercept = rep(NA,length(clstatlist)), 
                              slope = NA,
                              r2 = NA,
                              pvalue = NA)
bootvpdcoefhold <- bootmagcoefhold <- bootfreezecoefhold <- bootpptcoefhold
#
for (i in 1:length(clstatlist)) {
  # MAP
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxpptout <- outlierID(clhold[,c("resepx","ppt")])
  clhold <- clhold[which(resepxpptout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  # clhold <- clhold[which(clhold$resepxpptout == FALSE),] # removes the 
  # outliers so designated in the original combined data (uncomment to
  # use the original combined data; Note: these lines have been removed
  # in the code below; can probably remove this once everything is working)
  lmhold <- lm(resepx~ppt, clhold)
  bootpptcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,lmp(lmhold))
  #
  # Vapor pressure deficit
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxprismjjavpdout <- outlierID(clhold[,c("resepx","prismjjavpd")])
  clhold <- clhold[which(resepxprismjjavpdout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  lmhold <- lm(resepx~prismjjavpd, clhold)
  bootvpdcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,lmp(lmhold))
  #
  # Precipitation event magnitude
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmagnitudeout <- outlierID(clhold[,c("resepx","magnitude")])
  clhold <- clhold[which(resepxmagnitudeout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  lmhold <- lm(resepx~magnitude, clhold)
  bootmagcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,lmp(lmhold))
  #
  # Freezing frequency
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmintempfreezedaysout <- outlierID(clhold[,c("resepx",
                                                    "mintempfreezedays")])
  clhold <- clhold[which(resepxmintempfreezedaysout == FALSE),] # removes  
  # the outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  lmhold <- lm(resepx~mintempfreezedays, clhold)
  bootfreezecoefhold[i,] <- c(coef(lmhold),
                              summary(lmhold)$r.squared,lmp(lmhold))
}
#
# Save the clstatlist for the A horizons since this is changed for the 
# Ap and B below:
Aclstatlist <- clstatlist
#
Aclimatebootpptcoef600km <- bootpptcoefhold
Aclimatebootvpdcoef600km <- bootvpdcoefhold 
Aclimatebootmagcoef600km <- bootmagcoefhold 
Aclimatebootfreezecoef600km <- bootfreezecoefhold
#
# 
# Output the RESEPx results of the bootstrapped aggregated data separated by
# the appropriate minimum distance
Aclimatexagg600kmlist <- Aclstatlist # indices
Aclimatexagg600km1 <- myclimatehold[Aclstatlist[[1]],]
Aclimatexagg600km2 <- myclimatehold[Aclstatlist[[2]],]
Aclimatexagg600km3 <- myclimatehold[Aclstatlist[[3]],]
Aclimatexagg600km4 <- myclimatehold[Aclstatlist[[4]],]
Aclimatexagg600km5 <- myclimatehold[Aclstatlist[[5]],]
#
#
# Record the outlier information
#
# A horizons
#
# rep1 
x1 <- Aclimatexagg600km1$ppt
x2 <- Aclimatexagg600km1$prismjjavpd
x3 <- Aclimatexagg600km1$magnitude
x4 <- Aclimatexagg600km1$mintempfreezedays
resepx <- Aclimatexagg600km1$resepx
Aclimatexagg600km1out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Aclimatexagg600km1out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Aclimatexagg600km1out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Aclimatexagg600km1out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Aclimatexagg600km1out$resepxfreezeout <- outlierID(dhold)
#
# rep2
x1 <- Aclimatexagg600km2$ppt
x2 <- Aclimatexagg600km2$prismjjavpd
x3 <- Aclimatexagg600km2$magnitude
x4 <- Aclimatexagg600km2$mintempfreezedays
resepx <- Aclimatexagg600km2$resepx
Aclimatexagg600km2out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Aclimatexagg600km2out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Aclimatexagg600km2out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Aclimatexagg600km2out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Aclimatexagg600km2out$resepxfreezeout <- outlierID(dhold)
#
# rep3
x1 <- Aclimatexagg600km3$ppt
x2 <- Aclimatexagg600km3$prismjjavpd
x3 <- Aclimatexagg600km3$magnitude
x4 <- Aclimatexagg600km3$mintempfreezedays
resepx <- Aclimatexagg600km3$resepx
Aclimatexagg600km3out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Aclimatexagg600km3out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Aclimatexagg600km3out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Aclimatexagg600km3out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Aclimatexagg600km3out$resepxfreezeout <- outlierID(dhold)
#
# rep4
x1 <- Aclimatexagg600km4$ppt
x2 <- Aclimatexagg600km4$prismjjavpd
x3 <- Aclimatexagg600km4$magnitude
x4 <- Aclimatexagg600km4$mintempfreezedays
resepx <- Aclimatexagg600km4$resepx
Aclimatexagg600km4out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Aclimatexagg600km4out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Aclimatexagg600km4out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Aclimatexagg600km4out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Aclimatexagg600km4out$resepxfreezeout <- outlierID(dhold)
#
# rep5
x1 <- Aclimatexagg600km5$ppt
x2 <- Aclimatexagg600km5$prismjjavpd
x3 <- Aclimatexagg600km5$magnitude
x4 <- Aclimatexagg600km5$mintempfreezedays
resepx <- Aclimatexagg600km5$resepx
Aclimatexagg600km5out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Aclimatexagg600km5out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Aclimatexagg600km5out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Aclimatexagg600km5out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Aclimatexagg600km5out$resepxfreezeout <- outlierID(dhold)
#
#
setwd(oname)
write.table(data.frame(Aclimatexagg600km1[,1:24],
                       Aclimatexagg600km1out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatexagg300km1.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Aclimatexagg600km2[,1:24],
                       Aclimatexagg600km2out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatexagg300km2.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Aclimatexagg600km3[,1:24],
                       Aclimatexagg600km3out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatexagg300km3.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Aclimatexagg600km4[,1:24],
                       Aclimatexagg600km4out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatexagg300km4.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Aclimatexagg600km5[,1:24],
                       Aclimatexagg600km5out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatexagg300km5.csv',sep=""),
            sep=",",row.names=F)
setwd(wname)
#
#
#
# Subsetting the aggregated Ap horizon RESEPx data with a minimum distance
# ------------------------------------------------------------------------
#
# Switch the variables to generalize in the code below
mindistance <- 100 # km; for aggregated A horizons
mindistance <- mindistance*1000 # m
myclimatehold <- Apclimatex2
#
# Create a distance matrix to find the minimum distance in meters
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
distmat <- matrix(NA, nrow=nrow(cl),ncol=nrow(cl))
for (j in 1:nrow(cl)) { # columns in distmat
  for (i in 1:nrow(cl)) { # rows in distmat
    distmat[i,j] <- distHaversine(cl[j,],cl[i,]) # calculates the distance
    # on the ground in meters using the Haversine formula
  }
}
#
#
# Run a loop to randomly select climate stations that are separated by
# a minimum distance. The data for these stations are the aggregated A 
# horizon data after outliers were removed prior to Ruth's previous
# analyses. The loop randomly selects these climate stations by using
# different starting points and order of stations checked each time.
# The loop is run of N times following the code below:
#
N <- 1200 # recommended by Lander as a good number for a general
# bootstrap
clstatlist <- list(NULL) # will hold the results of the loop below
# (run N times) in each element; results are the distmat indices
# that are at the appropriate spacing
#
for (i in 1:N) {
  # Run a loop to select climate stations with a minimum distance
  rdistmatcolindex <- sample(1:ncol(distmat)) # randomizes the column index
  # for distmat so that it doesn't keep starting in the same place
  myindexhold <- sample(which(distmat[,rdistmatcolindex[1]] > mindistance)) 
  # forms the initial maximum subset of the data
  for (j in 2:length(myindexhold)) {
    mycheck <- which(distmat[,myindexhold[j]] < 
                       mindistance) # holds the distmat indices of the 
    # stations that are too close to the jth myindexhold station under 
    # consideration
    mycheck <- mycheck[-which(mycheck == myindexhold[j])] # removes the
    # distmat index of the jth myindexhold station under consideration
    myhold <- NULL # resets myhold
    #
    myhold <- ms(myindexhold,mycheck) # holds the indices of the 
    # myindexhold vector that should be removed
    #
    if (length(myhold) == 0) myhold <- NULL # no indices were found; 
    # reset myhold to NULL for the if statement below
    #
    if(!is.null(myhold)) myindexhold <- myindexhold[-myhold]
    if (j >= length(myindexhold)) break
  }
  clstatlist[[i]] <- myindexhold
}
#
#
# Check the distribution of all the indices to see how many of the climate
# stations are included in the analysis:
map('usa')
map("state", boundary = FALSE, col="gray", add = TRUE)
mycolors <- jet.colors(length(clstatlist))
for (i in 1:length(clstatlist)){
  points(lat~long,Apclimatex2[clstatlist[[i]],],cex=0.5,pch=21,
         bg=mycolors[i])
}
#
#
# Run the regressions for each of the appropriately-spaced subsets in
# clstatlist:
bootpptcoefhold <- data.frame(intercept = rep(NA,length(clstatlist)), 
                              slope = NA,
                              r2 = NA,
                              pvalue = NA)
bootvpdcoefhold <- bootmagcoefhold <- bootfreezecoefhold <- bootpptcoefhold
for (i in 1:length(clstatlist)) {
  # MAP
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxpptout <- outlierID(clhold[,c("resepx","ppt")])
  clhold <- clhold[which(resepxpptout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  lmhold <- lm(resepx~ppt, clhold)
  bootpptcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,lmp(lmhold))
  #
  # Vapor pressure deficit
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxprismjjavpdout <- outlierID(clhold[,c("resepx","prismjjavpd")])
  clhold <- clhold[which(resepxprismjjavpdout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  lmhold <- lm(resepx~prismjjavpd, clhold)
  bootvpdcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,lmp(lmhold))
  #
  # Precipitation event magnitude
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmagnitudeout <- outlierID(clhold[,c("resepx","magnitude")])
  clhold <- clhold[which(resepxmagnitudeout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  lmhold <- lm(resepx~magnitude, clhold)
  bootmagcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,lmp(lmhold))
  #
  # Freezing frequency
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmintempfreezedaysout <- outlierID(clhold[,c("resepx",
                                                    "mintempfreezedays")])
  clhold <- clhold[which(resepxmintempfreezedaysout == FALSE),] # removes  
  # the outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  lmhold <- lm(resepx~mintempfreezedays, clhold)
  bootfreezecoefhold[i,] <- c(coef(lmhold),
                              summary(lmhold)$r.squared,lmp(lmhold))
}
#
# Save the clstatlist for the Ap horizons since this is changed for the 
# A above and B below:
Apclstatlist <- clstatlist
#
Apclimatebootpptcoef220km <- bootpptcoefhold
Apclimatebootvpdcoef220km <- bootvpdcoefhold 
Apclimatebootmagcoef220km <- bootmagcoefhold 
Apclimatebootfreezecoef220km <- bootfreezecoefhold
#
# 
# Output the RESEPx results of the bootstrapped aggregated data separated by
# the appropriate minimum distance
Apclimatexagg220kmlist <- clstatlist # indices
Apclimatexagg220km1 <- myclimatehold[Apclstatlist[[1]],]
Apclimatexagg220km2 <- myclimatehold[Apclstatlist[[2]],]
Apclimatexagg220km3 <- myclimatehold[Apclstatlist[[3]],]
Apclimatexagg220km4 <- myclimatehold[Apclstatlist[[4]],]
Apclimatexagg220km5 <- myclimatehold[Apclstatlist[[5]],]
#
#
# Record the outliers:
#
# Ap horizons
#
# rep1 
x1 <- Apclimatexagg220km1$ppt
x2 <- Apclimatexagg220km1$prismjjavpd
x3 <- Apclimatexagg220km1$magnitude
x4 <- Apclimatexagg220km1$mintempfreezedays
resepx <- Apclimatexagg220km1$resepx
Apclimatexagg220km1out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                     resepxvpdout=FALSE,
                                     resepxmagout=FALSE,
                                     resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Apclimatexagg220km1out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Apclimatexagg220km1out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Apclimatexagg220km1out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Apclimatexagg220km1out$resepxfreezeout <- outlierID(dhold)
#
# rep2
x1 <- Apclimatexagg220km2$ppt
x2 <- Apclimatexagg220km2$prismjjavpd
x3 <- Apclimatexagg220km2$magnitude
x4 <- Apclimatexagg220km2$mintempfreezedays
resepx <- Apclimatexagg220km2$resepx
Apclimatexagg220km2out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                     resepxvpdout=FALSE,
                                     resepxmagout=FALSE,
                                     resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Apclimatexagg220km2out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Apclimatexagg220km2out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Apclimatexagg220km2out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Apclimatexagg220km2out$resepxfreezeout <- outlierID(dhold)
#
# rep3
x1 <- Apclimatexagg220km3$ppt
x2 <- Apclimatexagg220km3$prismjjavpd
x3 <- Apclimatexagg220km3$magnitude
x4 <- Apclimatexagg220km3$mintempfreezedays
resepx <- Apclimatexagg220km3$resepx
Apclimatexagg220km3out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                     resepxvpdout=FALSE,
                                     resepxmagout=FALSE,
                                     resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Apclimatexagg220km3out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Apclimatexagg220km3out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Apclimatexagg220km3out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Apclimatexagg220km3out$resepxfreezeout <- outlierID(dhold)
#
# rep4
x1 <- Apclimatexagg220km4$ppt
x2 <- Apclimatexagg220km4$prismjjavpd
x3 <- Apclimatexagg220km4$magnitude
x4 <- Apclimatexagg220km4$mintempfreezedays
resepx <- Apclimatexagg220km4$resepx
Apclimatexagg220km4out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                     resepxvpdout=FALSE,
                                     resepxmagout=FALSE,
                                     resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Apclimatexagg220km4out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Apclimatexagg220km4out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Apclimatexagg220km4out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Apclimatexagg220km4out$resepxfreezeout <- outlierID(dhold)
#
# rep5
x1 <- Apclimatexagg220km5$ppt
x2 <- Apclimatexagg220km5$prismjjavpd
x3 <- Apclimatexagg220km5$magnitude
x4 <- Apclimatexagg220km5$mintempfreezedays
resepx <- Apclimatexagg220km5$resepx
Apclimatexagg220km5out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                     resepxvpdout=FALSE,
                                     resepxmagout=FALSE,
                                     resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Apclimatexagg220km5out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Apclimatexagg220km5out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Apclimatexagg220km5out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Apclimatexagg220km5out$resepxfreezeout <- outlierID(dhold)
#
#
setwd(oname)
write.table(data.frame(Apclimatexagg220km1[,1:24],
                       Apclimatexagg220km1out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatexagg100km1.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Apclimatexagg220km2[,1:24],
                       Apclimatexagg220km2out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatexagg100km2.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Apclimatexagg220km3[,1:24],
                       Apclimatexagg220km3out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatexagg100km3.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Apclimatexagg220km4[,1:24],
                       Apclimatexagg220km4out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatexagg100km4.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Apclimatexagg220km5[,1:24],
                       Apclimatexagg220km5out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatexagg100km5.csv',sep=""),
            sep=",",row.names=F)
setwd(wname)
#
#
# Subsetting the aggregated B horizon RESEPx data with a minimum distance
# -----------------------------------------------------------------------
#
# Switch the variables to generalize in the code below
mindistance <- 420 # km; for aggregated A horizons
mindistance <- mindistance*1000 # m
myclimatehold <- Bclimatex2
#
# Create a distance matrix to find the minimum distance in meters
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
distmat <- matrix(NA, nrow=nrow(cl),ncol=nrow(cl))
for (j in 1:nrow(cl)) { # columns in distmat
  for (i in 1:nrow(cl)) { # rows in distmat
    distmat[i,j] <- distHaversine(cl[j,],cl[i,]) # calculates the distance
    # on the ground in meters using the Haversine formula
  }
}
#
#
# Run a loop to randomly select climate stations that are separated by
# a minimum distance. The data for these stations are the aggregated A 
# horizon data after outliers were removed prior to Ruth's previous
# analyses. The loop randomly selects these climate stations by using
# different starting points and order of stations checked each time.
# The loop is run of N times following the code below:
#
N <- 1200 # recommended by Lander as a good number for a general
# bootstrap
clstatlist <- list(NULL) # will hold the results of the loop below
# (run N times) in each element; results are the distmat indices
# that are at the appropriate spacing
#
for (i in 1:N) {
  # Run a loop to select climate stations with a minimum distance
  rdistmatcolindex <- sample(1:ncol(distmat)) # randomizes the column index
  # for distmat so that it doesn't keep starting in the same place
  myindexhold <- sample(which(distmat[,rdistmatcolindex[1]] > mindistance)) 
  # forms the initial maximum subset of the data
  for (j in 2:length(myindexhold)) {
    mycheck <- which(distmat[,myindexhold[j]] < 
                       mindistance) # holds the distmat indices of the 
    # stations that are too close to the jth myindexhold station under 
    # consideration
    mycheck <- mycheck[-which(mycheck == myindexhold[j])] # removes the
    # distmat index of the jth myindexhold station under consideration
    myhold <- NULL # resets myhold
    #
    myhold <- ms(myindexhold,mycheck) # holds the indices of the 
    # myindexhold vector that should be removed
    #
    if (length(myhold) == 0) myhold <- NULL # no indices were found; 
    # reset myhold to NULL for the if statement below
    #
    if(!is.null(myhold)) myindexhold <- myindexhold[-myhold]
    if (j >= length(myindexhold)) break
  }
  clstatlist[[i]] <- myindexhold
}
#
#
# Check the distribution of all the indices to see how many of the climate
# stations are included in the analysis:
map('usa')
map("state", boundary = FALSE, col="gray", add = TRUE)
mycolors <- jet.colors(length(clstatlist))
for (i in 1:length(clstatlist)){
  points(lat~long,Bclimatex2[clstatlist[[i]],],cex=0.5,pch=21,bg=mycolors[i])
}
#
#
# Run the regressions for each of the appropriately-spaced subsets in
# clstatlist:
bootpptcoefhold <- data.frame(intercept = rep(NA,length(clstatlist)), 
                              slope = NA,
                              r2 = NA,
                              pvalue = NA)
bootvpdcoefhold <- bootmagcoefhold <- bootfreezecoefhold <- bootpptcoefhold
for (i in 1:length(clstatlist)) {
  # MAP
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxpptout <- outlierID(clhold[,c("resepx","ppt")])
  clhold <- clhold[which(resepxpptout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  lmhold <- lm(resepx~ppt, clhold)
  bootpptcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,lmp(lmhold))
  #
  # Vapor pressure deficit
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxprismjjavpdout <- outlierID(clhold[,c("resepx","prismjjavpd")])
  clhold <- clhold[which(resepxprismjjavpdout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  lmhold <- lm(resepx~prismjjavpd, clhold)
  bootvpdcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,lmp(lmhold))
  #
  # Precipitation event magnitude
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmagnitudeout <- outlierID(clhold[,c("resepx","magnitude")])
  clhold <- clhold[which(resepxmagnitudeout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  lmhold <- lm(resepx~magnitude, clhold)
  bootmagcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,lmp(lmhold))
  #
  # Freezing frequency
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmintempfreezedaysout <- outlierID(clhold[,c("resepx",
                                                    "mintempfreezedays")])
  clhold <- clhold[which(resepxmintempfreezedaysout == FALSE),] # removes  
  # the outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  lmhold <- lm(resepx~mintempfreezedays, clhold)
  bootfreezecoefhold[i,] <- c(coef(lmhold),
                              summary(lmhold)$r.squared,lmp(lmhold))
}
#
# Save the clstatlist for the B horizons since this is changed for the 
# A and Ap above:
Bclstatlist <- clstatlist
#
Bclimatebootpptcoef350km <- bootpptcoefhold
Bclimatebootvpdcoef350km <- bootvpdcoefhold 
Bclimatebootmagcoef350km <- bootmagcoefhold 
Bclimatebootfreezecoef350km <- bootfreezecoefhold
#
#
# Output the RESEPx results of the bootstrapped aggregated data separated by
# the appropriate minimum distance
Bclimatexagg350kmlist <- Bclstatlist # indices
Bclimatexagg350km1 <- myclimatehold[Bclstatlist[[1]],]
Bclimatexagg350km2 <- myclimatehold[Bclstatlist[[2]],]
Bclimatexagg350km3 <- myclimatehold[Bclstatlist[[3]],]
Bclimatexagg350km4 <- myclimatehold[Bclstatlist[[4]],]
Bclimatexagg350km5 <- myclimatehold[Bclstatlist[[5]],]
#
#
# Record the outlier information
#
# B horizons
#
# rep1 
x1 <- Bclimatexagg350km1$ppt
x2 <- Bclimatexagg350km1$prismjjavpd
x3 <- Bclimatexagg350km1$magnitude
x4 <- Bclimatexagg350km1$mintempfreezedays
resepx <- Bclimatexagg350km1$resepx
Bclimatexagg350km1out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Bclimatexagg350km1out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Bclimatexagg350km1out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Bclimatexagg350km1out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Bclimatexagg350km1out$resepxfreezeout <- outlierID(dhold)
#
# rep2
x1 <- Bclimatexagg350km2$ppt
x2 <- Bclimatexagg350km2$prismjjavpd
x3 <- Bclimatexagg350km2$magnitude
x4 <- Bclimatexagg350km2$mintempfreezedays
resepx <- Bclimatexagg350km2$resepx
Bclimatexagg350km2out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Bclimatexagg350km2out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Bclimatexagg350km2out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Bclimatexagg350km2out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Bclimatexagg350km2out$resepxfreezeout <- outlierID(dhold)
#
# rep3
x1 <- Bclimatexagg350km3$ppt
x2 <- Bclimatexagg350km3$prismjjavpd
x3 <- Bclimatexagg350km3$magnitude
x4 <- Bclimatexagg350km3$mintempfreezedays
resepx <- Bclimatexagg350km3$resepx
Bclimatexagg350km3out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Bclimatexagg350km3out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Bclimatexagg350km3out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Bclimatexagg350km3out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Bclimatexagg350km3out$resepxfreezeout <- outlierID(dhold)
#
# rep4
x1 <- Bclimatexagg350km4$ppt
x2 <- Bclimatexagg350km4$prismjjavpd
x3 <- Bclimatexagg350km4$magnitude
x4 <- Bclimatexagg350km4$mintempfreezedays
resepx <- Bclimatexagg350km4$resepx
Bclimatexagg350km4out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Bclimatexagg350km4out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Bclimatexagg350km4out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Bclimatexagg350km4out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Bclimatexagg350km4out$resepxfreezeout <- outlierID(dhold)
#
# rep5
x1 <- Bclimatexagg350km5$ppt
x2 <- Bclimatexagg350km5$prismjjavpd
x3 <- Bclimatexagg350km5$magnitude
x4 <- Bclimatexagg350km5$mintempfreezedays
resepx <- Bclimatexagg350km5$resepx
Bclimatexagg350km5out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Bclimatexagg350km5out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Bclimatexagg350km5out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Bclimatexagg350km5out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Bclimatexagg350km5out$resepxfreezeout <- outlierID(dhold)
#
#
setwd(oname)
write.table(data.frame(Bclimatexagg350km1[,1:24],
                       Bclimatexagg350km1out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatexagg420km1.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Bclimatexagg350km2[,1:24],
                       Bclimatexagg350km2out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatexagg420km2.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Bclimatexagg350km3[,1:24],
                       Bclimatexagg350km3out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatexagg420km3.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Bclimatexagg350km4[,1:24],
                       Bclimatexagg350km4out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatexagg420km4.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Bclimatexagg350km5[,1:24],
                       Bclimatexagg350km5out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatexagg420km5.csv',sep=""),
            sep=",",row.names=F)
setwd(wname)
#
#
# Output the results of the bootstrapped aggregated data separated by
# the appropriate minimum distance
# -------------------------------------------------------------------
setwd(oname)
# Summary Data:
#
# A horizons
write.table(Aclimatebootpptcoef600km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatebootpptcoef600km.csv',sep=""),
            sep=",",row.names=F)
write.table(Aclimatebootvpdcoef600km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatebootvpdcoef600km.csv',sep=""),
            sep=",",row.names=F)
write.table(Aclimatebootmagcoef600km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatebootmagcoef600km.csv',sep=""),
            sep=",",row.names=F)
write.table(Aclimatebootfreezecoef600km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Aclimatebootfreezecoef600km.csv',sep=""),
            sep=",",row.names=F)
#
# Ap horizons
write.table(Apclimatebootpptcoef220km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatebootpptcoef220km.csv',sep=""),
            sep=",",row.names=F)
write.table(Apclimatebootvpdcoef220km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatebootvpdcoef220km.csv',sep=""),
            sep=",",row.names=F)
write.table(Apclimatebootmagcoef220km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatebootmagcoef220km.csv',sep=""),
            sep=",",row.names=F)
write.table(Apclimatebootfreezecoef220km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Apclimatebootfreezecoef220km.csv',sep=""),
            sep=",",row.names=F)
#
# B horizons
write.table(Bclimatebootpptcoef350km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatebootpptcoef350km.csv',sep=""),
            sep=",",row.names=F)
write.table(Bclimatebootvpdcoef350km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatebootvpdcoef350km.csv',sep=""),
            sep=",",row.names=F)
write.table(Bclimatebootmagcoef350km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatebootmagcoef350km.csv',sep=""),
            sep=",",row.names=F)
write.table(Bclimatebootfreezecoef350km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_Bclimatebootfreezecoef350km.csv',sep=""),
            sep=",",row.names=F)
setwd(wname)
#
#
#
# Aggregate the results of the regression coefficients
# ----------------------------------------------------
slopemeanmindist <- data.frame(horizon=c("A","Ap","B"),
                                  ppt=NA,
                                  vpd=NA,
                                  mag=NA,
                                  freeze=NA)
interceptmeanmindist <- slopemeanmindist
slopesdmindist <- slopemeanmindist
interceptsdmindist <- slopemeanmindist
slopePvaluemindist <- slopemeanmindist
interceptPvaluemindist <- slopemeanmindist
bootAlist <- list(Aclimatebootpptcoef600km,
                     Aclimatebootvpdcoef600km,
                     Aclimatebootmagcoef600km,
                     Aclimatebootfreezecoef600km)
bootAplist <- list(Apclimatebootpptcoef220km,
                   Apclimatebootvpdcoef220km,
                   Apclimatebootmagcoef220km,
                   Apclimatebootfreezecoef220km)
bootBlist <- list(Bclimatebootpptcoef350km,
                  Bclimatebootvpdcoef350km,
                  Bclimatebootmagcoef350km,
                  Bclimatebootfreezecoef350km)
for (i in 1:length(bootAlist)) {
  slopemeanmindist[1,i+1] <- mean(bootAlist[[i]][,"slope"])
  interceptmeanmindist[1,i+1] <- mean(bootAlist[[i]][,"intercept"])
  slopesdmindist[1,i+1] <- sd(bootAlist[[i]][,"slope"])
  interceptsdmindist[1,i+1] <- sd(bootAlist[[i]][,"intercept"])
  slopePvaluemindist[1,i+1] <- t.test(bootAlist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindist[1,i+1] <- t.test(bootAlist[[i]][,"intercept"],
                                          mu=0)$p.value
  slopemeanmindist[2,i+1] <- mean(bootAplist[[i]][,"slope"])
  interceptmeanmindist[2,i+1] <- mean(bootAplist[[i]][,"intercept"])
  slopesdmindist[2,i+1] <- sd(bootAplist[[i]][,"slope"])
  interceptsdmindist[2,i+1] <- sd(bootAplist[[i]][,"intercept"])
  slopePvaluemindist[2,i+1] <- t.test(bootAplist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindist[2,i+1] <- t.test(bootAplist[[i]][,"intercept"],
                                          mu=0)$p.value
  slopemeanmindist[3,i+1] <- mean(bootBlist[[i]][,"slope"])
  interceptmeanmindist[3,i+1] <- mean(bootBlist[[i]][,"intercept"])
  slopesdmindist[3,i+1] <- sd(bootBlist[[i]][,"slope"])
  interceptsdmindist[3,i+1] <- sd(bootBlist[[i]][,"intercept"])
  slopePvaluemindist[3,i+1] <- t.test(bootBlist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindist[3,i+1] <- t.test(bootBlist[[i]][,"intercept"],
                                          mu=0)$p.value
}
#
# Note: If this works, need to make plots of using these values and plot
# them to a new Figure 2 like figure.
# [insert code here]
#
#
# Output the summary of the results of the bootstrapped aggregated data 
# separated by the appropriate minimum distance
# ---------------------------------------------------------------------
setwd(oname)
write.table(slopemeanmindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_slopemeanmindist.csv',sep=""),
            sep=",",row.names=F)
write.table(interceptmeanmindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_interceptmeanmindist.csv',sep=""),
            sep=",",row.names=F)
write.table(slopesdmindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_slopesdmindist.csv',sep=""),
            sep=",",row.names=F)
write.table(interceptsdmindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_interceptsdmindist.csv',sep=""),
            sep=",",row.names=F)
write.table(slopePvaluemindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_slopePvaluemindist.csv',sep=""),
            sep=",",row.names=F)

write.table(interceptPvaluemindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_interceptPvaluemindist.csv',sep=""),
            sep=",",row.names=F)
setwd(wname)
#
# 
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Find the number of samples per station that is being used to create
# the aggregated data 5/7/18
# ---------------------------------------------------------------------------
Aclstatnum <- stats::aggregate(RESEPx ~ stnid, dAx, length)
names(Aclstatnum) <- c("stnid","N") # N is the number of resepx values for
# each station
Apclstatnum <- stats::aggregate(RESEPx ~ stnid, dApx, length)
names(Apclstatnum) <- c("stnid","N")
Bclstatnum <- stats::aggregate(RESEPx ~ stnid, dBx, length)
names(Bclstatnum) <- c("stnid","N")
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Rerun the bootstrapped estimates of the coefficients to account for 
# weighted regression (weights are the number of samples per station)
# 5/7/18
# ---------------------------------------------------------------------------
# From earlier but modified to account for these weights:
#
# Randomly sample the data so that samples are separated by some minimum 
# distance
# ----------------------------------------------------------------------
# Explanation: The minimum distances were set by Ruth's calculations in
# the previous analysis. Samples must be separated by at least 600 km
# for A horizons, 220 km for Ap horizons, and 350 km for B horizons. This
# is run using the already aggregated data and sampled many times to 
# estimate the parameters of the regression. A climate station with 
# aggregated data is chosen at random and the order of theothers are chosen 
# random in each iteration so that the minimum distancesfollow the results 
# in Ruth's analysis. Each realization of the subset with the appropriate 
# distance is used to estimate the regression parameters.
#
# What follows is large chunks of code for each horizon with the data
# outputted at the very end of the analysis.
#
# Note: there are several functions defined in the A horizon code below
# so those functions are not redefined in the Ap or B horizon code below
# that.
#
# Note: these have been updated (5/9/18) to 300, 100, and 420 km for
# A, Ap, and B horizons respectively; PLEASE NOTE!!!: Naming of the 
# data.frame follows the previous range values even though the ranges
# have been updated with Ruth's analysis on 5/9/18
#
# Subsetting the aggregated A horizon RESEPx data with a minimum distance
# -----------------------------------------------------------------------
#
# Switch the variables to generalize in the code below
mindistance <- 600 # km; for aggregated A horizons
mindistance <- mindistance*1000 # m
myclimatehold <- Aclimatex2
#
# Need the function distHaversine from the  package
require(geosphere)
#
# Create a distance matrix to find the minimum distance in meters
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
distmat <- matrix(NA, nrow=nrow(cl),ncol=nrow(cl))
for (j in 1:nrow(cl)) { # columns in distmat
  for (i in 1:nrow(cl)) { # rows in distmat
    distmat[i,j] <- distHaversine(cl[j,],cl[i,]) # calculates the distance
    # on the ground in meters using the Haversine formula
  }
}
#
#
# Function for finding the indices in a vector that correspond to multiple 
# conditions:
ms <- function(v,x,zequal="equal") {
  # v is a vector to search to see if the conditions are met
  # x is a vector of the conditions
  # zequal is either "equal", "greaterthan", or "lessthan"
  if (zequal == "equal") {
    out <- NULL
    for (i in 1:length(x)) {
      out <- c(out,which(v == x[i]))
    }
    out
  } else {
    if (zequal == "greaterthan") {
      out <- NULL
      for (i in 1:length(x)) {
        out <- c(out,which(v > x[i]))
      }
      out
    } else { # lessthan
      out <- NULL
      for (i in 1:length(x)) {
        out <- c(out,which(v < x[i]))
      }
      out
    }
  }
}
#
#
# Run a loop to randomly select climate stations that are separated by
# a minimum distance. The data for these stations are the aggregated A 
# horizon data after outliers were removed prior to Ruth's previous
# analyses. The loop randomly selects these climate stations by using
# different starting points and order of stations checked each time.
# The loop is run of N times following the code below:
#
N <- 8000 # 1200 is recommended by Lander as a good number for a general
# bootstrap; because the code (farther below) tests the spatial structure
# of the response variable in each resampling, here I use a larger number
# to get at least 1200 with all the NAs
clstatlist <- list(NULL) # will hold the results of the loop below
# (run N times) in each element; results are the distmat indices
# that are at the appropriate spacing
#
for (i in 1:N) {
  # Run a loop to select climate stations with a minimum distance
  rdistmatcolindex <- sample(1:ncol(distmat)) # randomizes the column index
  # for distmat so that it doesn't keep starting in the same place
  myindexhold <- sample(which(distmat[,rdistmatcolindex[1]] > mindistance))
  # forms the initial maximum subset of the data
  for (j in 2:length(myindexhold)) {
    mycheck <- which(distmat[,myindexhold[j]] < 
                       mindistance) # holds the distmat indices of the 
    # stations that are too close to the jth myindexhold station under 
    # consideration
    mycheck <- mycheck[-which(mycheck == myindexhold[j])] # removes the
    # distmat index of the jth myindexhold station under consideration
    myhold <- NULL # resets myhold
    #
    myhold <- ms(myindexhold,mycheck) # holds the indices of the 
    # myindexhold vector that should be removed
    #
    if (length(myhold) == 0) myhold <- NULL # no indices were found; 
    # reset myhold to NULL for the if statement below
    #
    if(!is.null(myhold)) myindexhold <- myindexhold[-myhold]
    if (j >= length(myindexhold)) break
  }
  clstatlist[[i]] <- myindexhold
}
#
#
# Check the distribution of all the indices to see how many of the climate
# stations are included in the analysis (bootstrapped so should be most of them):
map('usa')
map("state", boundary = FALSE, col="gray", add = TRUE)
mycolors <- jet.colors(length(clstatlist))
for (i in 1:length(clstatlist)){
  points(lat~long,Aclimatex2[clstatlist[[i]],],cex=0.5,pch=21,bg=mycolors[i])
}
#
#
# Function to extract a pvalue from an lm object
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
#
#
# Create a function to check for spatial structure in the response
checkMoransI <- function(df,response=names(df)[ncol(df)],x="long",
                         y="lat",alpha=0.05,neighb="idw",k=5) {
  require(ape)
  # Follows the example in: 
  # https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/
  # for idw
  #
  # df = data.frame
  # response = the df column name that spatial structure is being 
  # evaluated; default is the last column name in df; need to provide
  # a character vector of length 1; will also take an actual numeric 
  # vector representing the values for which spatial structure is being
  # evaluated
  # x,y = df column names that provide the geographic coordinates of 
  # each of the rows in df; need to provide strings if not "long" and 
  # "lat"
  # alpha = alpha level to decide signifance in the Moran's I test
  # neighb = type of neighborhood weights matrix; either "knn" for 
  # k nearest neighbor (k must be supplied) or "idw" for inverse
  # distance weight
  # k number of nearest neighbors; ignored if neighb = "idw"
  #
  if (neighb=="idw") {
    df.dists <- as.matrix(dist(cbind(df[,x], df[,y])))
    df.dists.inv <- 1/df.dists
    diag(df.dists.inv) <- 0
    W <- df.dists.inv # weights matrix
  } else { # "knn"
    require(sp); require(spdep)
    df <- as.data.frame(df)
    df.sp <- SpatialPointsDataFrame(c(df[,c(x,y)]),data = df)
    coords <- coordinates(df.sp)
    df.knn <- knearneigh(coords,k=k)
    df.nb <- knn2nb(df.knn)
    W <- nb2mat(df.nb) # weights matrix
    attributes(W)$call <- NULL
  }
  #
  out <- list(check=NULL,I=NULL,Ipvalue=NULL)
  if (is.character(response)) { # name of the response variable column
    moransI <- Moran.I(df[,response], W)
  } else { # actual values are provided
    moransI <- Moran.I(response, W)
  }
  #
  if(moransI$p.value > alpha) {
    check <- "No spatial structure"
  } else {
    check <- "Likely spatial structure"
  }
  out$check <- check
  out$I <- moransI$observed
  out$Ipvalue <- moransI$p.value
  out
}
#
#
# Run the regressions for each of the appropriately-spaced subsets in
# clstatlist:
bootpptcoefhold <- data.frame(intercept = rep(NA,length(clstatlist)), 
                              slope = NA,
                              r2 = NA,
                              pvalue = NA,
                              I = NA,
                              Ipvalue = NA)
bootvpdcoefhold <- bootmagcoefhold <- bootfreezecoefhold <- bootpptcoefhold
#
for (i in 1:length(clstatlist)) {
  # MAP
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxpptout <- outlierID(clhold[,c("resepx","ppt")])
  clhold <- clhold[which(resepxpptout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Aclstatnum,by="stnid")
  lmhold <- lm(resepx~ppt, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootpptcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                             lmp(lmhold),moransI$I,moransI$Ipvalue)
  #
  # Vapor pressure deficit
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxprismjjavpdout <- outlierID(clhold[,c("resepx","prismjjavpd")])
  clhold <- clhold[which(resepxprismjjavpdout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Aclstatnum,by="stnid")
  lmhold <- lm(resepx~prismjjavpd, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootvpdcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                           lmp(lmhold),moransI$I,moransI$Ipvalue)
  #
  # Precipitation event magnitude
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmagnitudeout <- outlierID(clhold[,c("resepx","magnitude")])
  clhold <- clhold[which(resepxmagnitudeout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Aclstatnum,by="stnid")
  lmhold <- lm(resepx~magnitude, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootmagcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                           lmp(lmhold),moransI$I,moransI$Ipvalue)
  #
  # Freezing frequency
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmintempfreezedaysout <- outlierID(clhold[,c("resepx",
                                                    "mintempfreezedays")])
  clhold <- clhold[which(resepxmintempfreezedaysout == FALSE),] # removes  
  # the outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Aclstatnum,by="stnid")
  lmhold <- lm(resepx~mintempfreezedays, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootfreezecoefhold[i,] <- c(coef(lmhold),
                              summary(lmhold)$r.squared,lmp(lmhold),
                              moransI$I,moransI$Ipvalue)
}
#
# Save the clstatlist for the A horizons since this is changed for the 
# Ap and B below:
Aclstatlist <- clstatlist
#
Aclimatebootpptcoef600km <- bootpptcoefhold
Aclimatebootvpdcoef600km <- bootvpdcoefhold 
Aclimatebootmagcoef600km <- bootmagcoefhold 
Aclimatebootfreezecoef600km <- bootfreezecoefhold
#
# 
# Output the RESEPx results of the bootstrapped aggregated data separated by
# the appropriate minimum distance
Aclimatexagg600kmlist <- Aclstatlist # indices
Aclimatexagg600km1 <- myclimatehold[Aclstatlist[[1]],]
Aclimatexagg600km2 <- myclimatehold[Aclstatlist[[2]],]
Aclimatexagg600km3 <- myclimatehold[Aclstatlist[[3]],]
Aclimatexagg600km4 <- myclimatehold[Aclstatlist[[4]],]
Aclimatexagg600km5 <- myclimatehold[Aclstatlist[[5]],]
#
#
# Record the outlier information
#
# A horizons
#
# rep1 
x1 <- Aclimatexagg600km1$ppt
x2 <- Aclimatexagg600km1$prismjjavpd
x3 <- Aclimatexagg600km1$magnitude
x4 <- Aclimatexagg600km1$mintempfreezedays
resepx <- Aclimatexagg600km1$resepx
Aclimatexagg600km1out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Aclimatexagg600km1out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Aclimatexagg600km1out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Aclimatexagg600km1out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Aclimatexagg600km1out$resepxfreezeout <- outlierID(dhold)
#
# rep2
x1 <- Aclimatexagg600km2$ppt
x2 <- Aclimatexagg600km2$prismjjavpd
x3 <- Aclimatexagg600km2$magnitude
x4 <- Aclimatexagg600km2$mintempfreezedays
resepx <- Aclimatexagg600km2$resepx
Aclimatexagg600km2out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Aclimatexagg600km2out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Aclimatexagg600km2out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Aclimatexagg600km2out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Aclimatexagg600km2out$resepxfreezeout <- outlierID(dhold)
#
# rep3
x1 <- Aclimatexagg600km3$ppt
x2 <- Aclimatexagg600km3$prismjjavpd
x3 <- Aclimatexagg600km3$magnitude
x4 <- Aclimatexagg600km3$mintempfreezedays
resepx <- Aclimatexagg600km3$resepx
Aclimatexagg600km3out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Aclimatexagg600km3out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Aclimatexagg600km3out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Aclimatexagg600km3out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Aclimatexagg600km3out$resepxfreezeout <- outlierID(dhold)
#
# rep4
x1 <- Aclimatexagg600km4$ppt
x2 <- Aclimatexagg600km4$prismjjavpd
x3 <- Aclimatexagg600km4$magnitude
x4 <- Aclimatexagg600km4$mintempfreezedays
resepx <- Aclimatexagg600km4$resepx
Aclimatexagg600km4out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Aclimatexagg600km4out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Aclimatexagg600km4out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Aclimatexagg600km4out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Aclimatexagg600km4out$resepxfreezeout <- outlierID(dhold)
#
# rep5
x1 <- Aclimatexagg600km5$ppt
x2 <- Aclimatexagg600km5$prismjjavpd
x3 <- Aclimatexagg600km5$magnitude
x4 <- Aclimatexagg600km5$mintempfreezedays
resepx <- Aclimatexagg600km5$resepx
Aclimatexagg600km5out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Aclimatexagg600km5out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Aclimatexagg600km5out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Aclimatexagg600km5out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Aclimatexagg600km5out$resepxfreezeout <- outlierID(dhold)
#
#
setwd(oname)
# w = weighted regression
write.table(data.frame(Aclimatexagg600km1[,1:24],
                       Aclimatexagg600km1out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wAclimatexagg600km1.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Aclimatexagg600km2[,1:24],
                       Aclimatexagg600km2out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wAclimatexagg600km2.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Aclimatexagg600km3[,1:24],
                       Aclimatexagg600km3out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wAclimatexagg600km3.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Aclimatexagg600km4[,1:24],
                       Aclimatexagg600km4out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wAclimatexagg600km4.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Aclimatexagg600km5[,1:24],
                       Aclimatexagg600km5out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wAclimatexagg600km5.csv',sep=""),
            sep=",",row.names=F)
setwd(wname)
#
#
#
# Subsetting the aggregated Ap horizon RESEPx data with a minimum distance
# ------------------------------------------------------------------------
#
# Switch the variables to generalize in the code below
mindistance <- 200 # km; for aggregated Ap horizons
mindistance <- mindistance*1000 # m
myclimatehold <- Apclimatex2
#
# Create a distance matrix to find the minimum distance in meters
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
distmat <- matrix(NA, nrow=nrow(cl),ncol=nrow(cl))
for (j in 1:nrow(cl)) { # columns in distmat
  for (i in 1:nrow(cl)) { # rows in distmat
    distmat[i,j] <- distHaversine(cl[j,],cl[i,]) # calculates the distance
    # on the ground in meters using the Haversine formula
  }
}
#
#
# Run a loop to randomly select climate stations that are separated by
# a minimum distance. The data for these stations are the aggregated A 
# horizon data after outliers were removed prior to Ruth's previous
# analyses. The loop randomly selects these climate stations by using
# different starting points and order of stations checked each time.
# The loop is run of N times following the code below:
#
clstatlist <- list(NULL) # will hold the results of the loop below
# (run N times) in each element; results are the distmat indices
# that are at the appropriate spacing
#
for (i in 1:N) {
  # Run a loop to select climate stations with a minimum distance
  rdistmatcolindex <- sample(1:ncol(distmat)) # randomizes the column index
  # for distmat so that it doesn't keep starting in the same place
  myindexhold <- sample(which(distmat[,rdistmatcolindex[1]] > mindistance)) 
  # forms the initial maximum subset of the data
  for (j in 2:length(myindexhold)) {
    mycheck <- which(distmat[,myindexhold[j]] < 
                       mindistance) # holds the distmat indices of the 
    # stations that are too close to the jth myindexhold station under 
    # consideration
    mycheck <- mycheck[-which(mycheck == myindexhold[j])] # removes the
    # distmat index of the jth myindexhold station under consideration
    myhold <- NULL # resets myhold
    #
    myhold <- ms(myindexhold,mycheck) # holds the indices of the 
    # myindexhold vector that should be removed
    #
    if (length(myhold) == 0) myhold <- NULL # no indices were found; 
    # reset myhold to NULL for the if statement below
    #
    if(!is.null(myhold)) myindexhold <- myindexhold[-myhold]
    if (j >= length(myindexhold)) break
  }
  clstatlist[[i]] <- myindexhold
}
#
#
# Check the distribution of all the indices to see how many of the climate
# stations are included in the analysis:
map('usa')
map("state", boundary = FALSE, col="gray", add = TRUE)
mycolors <- jet.colors(length(clstatlist))
for (i in 1:length(clstatlist)){
  points(lat~long,Aclimatex2[clstatlist[[i]],],cex=0.5,pch=21,bg=mycolors[i])
}
#
#
# Run the regressions for each of the appropriately-spaced subsets in
# clstatlist:
bootpptcoefhold <- data.frame(intercept = rep(NA,length(clstatlist)), 
                              slope = NA,
                              r2 = NA,
                              pvalue = NA,
                              I = NA,
                              Ipvalue = NA)
bootvpdcoefhold <- bootmagcoefhold <- bootfreezecoefhold <- bootpptcoefhold
for (i in 1:length(clstatlist)) {
  # MAP
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxpptout <- outlierID(clhold[,c("resepx","ppt")])
  clhold <- clhold[which(resepxpptout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Apclstatnum,by="stnid")
  lmhold <- lm(resepx~ppt, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootpptcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                           lmp(lmhold),moransI$I,moransI$Ipvalue)
  #
  # Vapor pressure deficit
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxprismjjavpdout <- outlierID(clhold[,c("resepx","prismjjavpd")])
  clhold <- clhold[which(resepxprismjjavpdout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Apclstatnum,by="stnid")
  lmhold <- lm(resepx~prismjjavpd, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootvpdcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                           lmp(lmhold),moransI$I,moransI$Ipvalue)
  #
  # Precipitation event magnitude
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmagnitudeout <- outlierID(clhold[,c("resepx","magnitude")])
  clhold <- clhold[which(resepxmagnitudeout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Apclstatnum,by="stnid")
  lmhold <- lm(resepx~magnitude, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootmagcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                           lmp(lmhold),moransI$I,moransI$Ipvalue)
  #
  # Freezing frequency
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmintempfreezedaysout <- outlierID(clhold[,c("resepx",
                                                    "mintempfreezedays")])
  clhold <- clhold[which(resepxmintempfreezedaysout == FALSE),] # removes  
  # the outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Apclstatnum,by="stnid")
  lmhold <- lm(resepx~mintempfreezedays, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootfreezecoefhold[i,] <- c(coef(lmhold),
                              summary(lmhold)$r.squared,
                              lmp(lmhold),moransI$I,moransI$Ipvalue)
}
#
# Save the clstatlist for the Ap horizons since this is changed for the 
# A above and B below:
Apclstatlist <- clstatlist
#
Apclimatebootpptcoef220km <- bootpptcoefhold
Apclimatebootvpdcoef220km <- bootvpdcoefhold 
Apclimatebootmagcoef220km <- bootmagcoefhold 
Apclimatebootfreezecoef220km <- bootfreezecoefhold
#
# 
# Output the RESEPx results of the bootstrapped aggregated data separated by
# the appropriate minimum distance
Apclimatexagg220kmlist <- clstatlist # indices
Apclimatexagg220km1 <- myclimatehold[Apclstatlist[[1]],]
Apclimatexagg220km2 <- myclimatehold[Apclstatlist[[2]],]
Apclimatexagg220km3 <- myclimatehold[Apclstatlist[[3]],]
Apclimatexagg220km4 <- myclimatehold[Apclstatlist[[4]],]
Apclimatexagg220km5 <- myclimatehold[Apclstatlist[[5]],]
#
#
# Record the outliers:
#
# Ap horizons
#
# rep1 
x1 <- Apclimatexagg220km1$ppt
x2 <- Apclimatexagg220km1$prismjjavpd
x3 <- Apclimatexagg220km1$magnitude
x4 <- Apclimatexagg220km1$mintempfreezedays
resepx <- Apclimatexagg220km1$resepx
Apclimatexagg220km1out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                     resepxvpdout=FALSE,
                                     resepxmagout=FALSE,
                                     resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Apclimatexagg220km1out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Apclimatexagg220km1out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Apclimatexagg220km1out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Apclimatexagg220km1out$resepxfreezeout <- outlierID(dhold)
#
# rep2
x1 <- Apclimatexagg220km2$ppt
x2 <- Apclimatexagg220km2$prismjjavpd
x3 <- Apclimatexagg220km2$magnitude
x4 <- Apclimatexagg220km2$mintempfreezedays
resepx <- Apclimatexagg220km2$resepx
Apclimatexagg220km2out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                     resepxvpdout=FALSE,
                                     resepxmagout=FALSE,
                                     resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Apclimatexagg220km2out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Apclimatexagg220km2out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Apclimatexagg220km2out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Apclimatexagg220km2out$resepxfreezeout <- outlierID(dhold)
#
# rep3
x1 <- Apclimatexagg220km3$ppt
x2 <- Apclimatexagg220km3$prismjjavpd
x3 <- Apclimatexagg220km3$magnitude
x4 <- Apclimatexagg220km3$mintempfreezedays
resepx <- Apclimatexagg220km3$resepx
Apclimatexagg220km3out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                     resepxvpdout=FALSE,
                                     resepxmagout=FALSE,
                                     resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Apclimatexagg220km3out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Apclimatexagg220km3out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Apclimatexagg220km3out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Apclimatexagg220km3out$resepxfreezeout <- outlierID(dhold)
#
# rep4
x1 <- Apclimatexagg220km4$ppt
x2 <- Apclimatexagg220km4$prismjjavpd
x3 <- Apclimatexagg220km4$magnitude
x4 <- Apclimatexagg220km4$mintempfreezedays
resepx <- Apclimatexagg220km4$resepx
Apclimatexagg220km4out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                     resepxvpdout=FALSE,
                                     resepxmagout=FALSE,
                                     resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Apclimatexagg220km4out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Apclimatexagg220km4out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Apclimatexagg220km4out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Apclimatexagg220km4out$resepxfreezeout <- outlierID(dhold)
#
# rep5
x1 <- Apclimatexagg220km5$ppt
x2 <- Apclimatexagg220km5$prismjjavpd
x3 <- Apclimatexagg220km5$magnitude
x4 <- Apclimatexagg220km5$mintempfreezedays
resepx <- Apclimatexagg220km5$resepx
Apclimatexagg220km5out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                     resepxvpdout=FALSE,
                                     resepxmagout=FALSE,
                                     resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Apclimatexagg220km5out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Apclimatexagg220km5out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Apclimatexagg220km5out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Apclimatexagg220km5out$resepxfreezeout <- outlierID(dhold)
#
#
setwd(oname)
# w = weighted regression
write.table(data.frame(Apclimatexagg220km1[,1:24],
                       Apclimatexagg220km1out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wApclimatexagg200km1.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Apclimatexagg220km2[,1:24],
                       Apclimatexagg220km2out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wApclimatexagg200km2.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Apclimatexagg220km3[,1:24],
                       Apclimatexagg220km3out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wApclimatexagg200km3.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Apclimatexagg220km4[,1:24],
                       Apclimatexagg220km4out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wApclimatexagg200km4.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Apclimatexagg220km5[,1:24],
                       Apclimatexagg220km5out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wApclimatexagg200km5.csv',sep=""),
            sep=",",row.names=F)
setwd(wname)
#
#
# Subsetting the aggregated B horizon RESEPx data with a minimum distance
# -----------------------------------------------------------------------
#
# Switch the variables to generalize in the code below
mindistance <- 840 # km; for aggregated B horizons
mindistance <- mindistance*1000 # m
myclimatehold <- Bclimatex2
#
# Create a distance matrix to find the minimum distance in meters
cl <- as.matrix(data.frame(lon=myclimatehold$long,lat=myclimatehold$lat))
distmat <- matrix(NA, nrow=nrow(cl),ncol=nrow(cl))
for (j in 1:nrow(cl)) { # columns in distmat
  for (i in 1:nrow(cl)) { # rows in distmat
    distmat[i,j] <- distHaversine(cl[j,],cl[i,]) # calculates the distance
    # on the ground in meters using the Haversine formula
  }
}
#
#
# Run a loop to randomly select climate stations that are separated by
# a minimum distance. The data for these stations are the aggregated A 
# horizon data after outliers were removed prior to Ruth's previous
# analyses. The loop randomly selects these climate stations by using
# different starting points and order of stations checked each time.
# The loop is run of N times following the code below:
#
clstatlist <- list(NULL) # will hold the results of the loop below
# (run N times) in each element; results are the distmat indices
# that are at the appropriate spacing
#
for (i in 1:N) {
  # Run a loop to select climate stations with a minimum distance
  rdistmatcolindex <- sample(1:ncol(distmat)) # randomizes the column index
  # for distmat so that it doesn't keep starting in the same place
  myindexhold <- sample(which(distmat[,rdistmatcolindex[1]] > mindistance)) 
  # forms the initial maximum subset of the data
  for (j in 2:length(myindexhold)) {
    mycheck <- which(distmat[,myindexhold[j]] < 
                       mindistance) # holds the distmat indices of the 
    # stations that are too close to the jth myindexhold station under 
    # consideration
    mycheck <- mycheck[-which(mycheck == myindexhold[j])] # removes the
    # distmat index of the jth myindexhold station under consideration
    myhold <- NULL # resets myhold
    #
    myhold <- ms(myindexhold,mycheck) # holds the indices of the 
    # myindexhold vector that should be removed
    #
    if (length(myhold) == 0) myhold <- NULL # no indices were found; 
    # reset myhold to NULL for the if statement below
    #
    if(!is.null(myhold)) myindexhold <- myindexhold[-myhold]
    if (j >= length(myindexhold)) break
  }
  clstatlist[[i]] <- myindexhold
}
#
#
# Check the distribution of all the indices to see how many of the climate
# stations are included in the analysis:
map('usa')
map("state", boundary = FALSE, col="gray", add = TRUE)
mycolors <- jet.colors(length(clstatlist))
for (i in 1:length(clstatlist)){
  points(lat~long,Aclimatex2[clstatlist[[i]],],cex=0.5,pch=21,bg=mycolors[i])
}
#
#
# Run the regressions for each of the appropriately-spaced subsets in
# clstatlist:
bootpptcoefhold <- data.frame(intercept = rep(NA,length(clstatlist)), 
                              slope = NA,
                              r2 = NA,
                              pvalue = NA,
                              I = NA,
                              Ipvalue = NA)
bootvpdcoefhold <- bootmagcoefhold <- bootfreezecoefhold <- bootpptcoefhold
for (i in 1:length(clstatlist)) {
  # MAP
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxpptout <- outlierID(clhold[,c("resepx","ppt")])
  clhold <- clhold[which(resepxpptout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Bclstatnum,by="stnid")
  lmhold <- lm(resepx~ppt, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootpptcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                           lmp(lmhold),moransI$I,moransI$Ipvalue)
  #
  # Vapor pressure deficit
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxprismjjavpdout <- outlierID(clhold[,c("resepx","prismjjavpd")])
  clhold <- clhold[which(resepxprismjjavpdout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Bclstatnum,by="stnid")
  lmhold <- lm(resepx~prismjjavpd, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootvpdcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                           lmp(lmhold),moransI$I,moransI$Ipvalue)
  #
  # Precipitation event magnitude
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmagnitudeout <- outlierID(clhold[,c("resepx","magnitude")])
  clhold <- clhold[which(resepxmagnitudeout == FALSE),] # removes the 
  # outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Bclstatnum,by="stnid")
  lmhold <- lm(resepx~magnitude, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootmagcoefhold[i,] <- c(coef(lmhold),summary(lmhold)$r.squared,
                           lmp(lmhold),moransI$I,moransI$Ipvalue)
  #
  # Freezing frequency
  clhold <- myclimatehold[clstatlist[[i]],]
  resepxmintempfreezedaysout <- outlierID(clhold[,c("resepx",
                                                    "mintempfreezedays")])
  clhold <- clhold[which(resepxmintempfreezedaysout == FALSE),] # removes  
  # the outliers so designated by the outlierID function; not saving this
  # information; using it only for these regressions; any checks of these
  # data will need to rerun the outlierID code
  clhold <- merge(clhold,Bclstatnum,by="stnid")
  lmhold <- lm(resepx~mintempfreezedays, clhold, weights = clhold$N)
  moransI <- checkMoransI(clhold,response=residuals(lmhold))
  bootfreezecoefhold[i,] <- c(coef(lmhold),
                              summary(lmhold)$r.squared,
                              lmp(lmhold),moransI$I,moransI$Ipvalue)
}
#
# Save the clstatlist for the B horizons since this is changed for the 
# A and Ap above:
Bclstatlist <- clstatlist
#
Bclimatebootpptcoef350km <- bootpptcoefhold
Bclimatebootvpdcoef350km <- bootvpdcoefhold 
Bclimatebootmagcoef350km <- bootmagcoefhold 
Bclimatebootfreezecoef350km <- bootfreezecoefhold
#
#
# Output the RESEPx results of the bootstrapped aggregated data separated by
# the appropriate minimum distance
Bclimatexagg350kmlist <- Bclstatlist # indices
Bclimatexagg350km1 <- myclimatehold[Bclstatlist[[1]],]
Bclimatexagg350km2 <- myclimatehold[Bclstatlist[[2]],]
Bclimatexagg350km3 <- myclimatehold[Bclstatlist[[3]],]
Bclimatexagg350km4 <- myclimatehold[Bclstatlist[[4]],]
Bclimatexagg350km5 <- myclimatehold[Bclstatlist[[5]],]
#
#
# Record the outlier information
#
# B horizons
#
# rep1 
x1 <- Bclimatexagg350km1$ppt
x2 <- Bclimatexagg350km1$prismjjavpd
x3 <- Bclimatexagg350km1$magnitude
x4 <- Bclimatexagg350km1$mintempfreezedays
resepx <- Bclimatexagg350km1$resepx
Bclimatexagg350km1out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Bclimatexagg350km1out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Bclimatexagg350km1out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Bclimatexagg350km1out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Bclimatexagg350km1out$resepxfreezeout <- outlierID(dhold)
#
# rep2
x1 <- Bclimatexagg350km2$ppt
x2 <- Bclimatexagg350km2$prismjjavpd
x3 <- Bclimatexagg350km2$magnitude
x4 <- Bclimatexagg350km2$mintempfreezedays
resepx <- Bclimatexagg350km2$resepx
Bclimatexagg350km2out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Bclimatexagg350km2out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Bclimatexagg350km2out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Bclimatexagg350km2out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Bclimatexagg350km2out$resepxfreezeout <- outlierID(dhold)
#
# rep3
x1 <- Bclimatexagg350km3$ppt
x2 <- Bclimatexagg350km3$prismjjavpd
x3 <- Bclimatexagg350km3$magnitude
x4 <- Bclimatexagg350km3$mintempfreezedays
resepx <- Bclimatexagg350km3$resepx
Bclimatexagg350km3out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Bclimatexagg350km3out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Bclimatexagg350km3out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Bclimatexagg350km3out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Bclimatexagg350km3out$resepxfreezeout <- outlierID(dhold)
#
# rep4
x1 <- Bclimatexagg350km4$ppt
x2 <- Bclimatexagg350km4$prismjjavpd
x3 <- Bclimatexagg350km4$magnitude
x4 <- Bclimatexagg350km4$mintempfreezedays
resepx <- Bclimatexagg350km4$resepx
Bclimatexagg350km4out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Bclimatexagg350km4out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Bclimatexagg350km4out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Bclimatexagg350km4out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Bclimatexagg350km4out$resepxfreezeout <- outlierID(dhold)
#
# rep5
x1 <- Bclimatexagg350km5$ppt
x2 <- Bclimatexagg350km5$prismjjavpd
x3 <- Bclimatexagg350km5$magnitude
x4 <- Bclimatexagg350km5$mintempfreezedays
resepx <- Bclimatexagg350km5$resepx
Bclimatexagg350km5out <- data.frame(resepxpptout=rep(FALSE,length(x1)),
                                    resepxvpdout=FALSE,
                                    resepxmagout=FALSE,
                                    resepxfreezeout=FALSE)
dhold <- data.frame(x=x1,resepx=resepx)
Bclimatexagg350km5out$resepxpptout <- outlierID(dhold)
dhold <- data.frame(x=x2,resepx=resepx)
Bclimatexagg350km5out$resepxvpdout <- outlierID(dhold)
dhold <- data.frame(x=x3,resepx=resepx)
Bclimatexagg350km5out$resepxmagout <- outlierID(dhold)
dhold <- data.frame(x=x4,resepx=resepx)
Bclimatexagg350km5out$resepxfreezeout <- outlierID(dhold)
#
#
setwd(oname)
# w = weighted regression
write.table(data.frame(Bclimatexagg350km1[,1:24],
                       Bclimatexagg350km1out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wBclimatexagg840km1.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Bclimatexagg350km2[,1:24],
                       Bclimatexagg350km2out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wBclimatexagg840km2.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Bclimatexagg350km3[,1:24],
                       Bclimatexagg350km3out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wBclimatexagg840km3.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Bclimatexagg350km4[,1:24],
                       Bclimatexagg350km4out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wBclimatexagg840km4.csv',sep=""),
            sep=",",row.names=F)
write.table(data.frame(Bclimatexagg350km5[,1:24],
                       Bclimatexagg350km5out),
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wBclimatexagg840km5.csv',sep=""),
            sep=",",row.names=F)
setwd(wname)
#
#
# Output the results of the bootstrapped aggregated data separated by
# the appropriate minimum distance
# -------------------------------------------------------------------
setwd(oname)
# w = weighted regression
# Summary Data:
#
# A horizons
write.table(Aclimatebootpptcoef600km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wAclimatebootpptcoef600km.csv',sep=""),
            sep=",",row.names=F)
write.table(Aclimatebootvpdcoef600km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wAclimatebootvpdcoef600km.csv',sep=""),
            sep=",",row.names=F)
write.table(Aclimatebootmagcoef600km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wAclimatebootmagcoef600km.csv',sep=""),
            sep=",",row.names=F)
write.table(Aclimatebootfreezecoef600km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wAclimatebootfreezecoef600km.csv',sep=""),
            sep=",",row.names=F)
#
# Ap horizons
write.table(Apclimatebootpptcoef220km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wApclimatebootpptcoef200km.csv',sep=""),
            sep=",",row.names=F)
write.table(Apclimatebootvpdcoef220km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wApclimatebootvpdcoef200km.csv',sep=""),
            sep=",",row.names=F)
write.table(Apclimatebootmagcoef220km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wApclimatebootmagcoef200km.csv',sep=""),
            sep=",",row.names=F)
write.table(Apclimatebootfreezecoef220km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wApclimatebootfreezecoef200km.csv',sep=""),
            sep=",",row.names=F)
#
# B horizons
write.table(Bclimatebootpptcoef350km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wBclimatebootpptcoef840km.csv',sep=""),
            sep=",",row.names=F)
write.table(Bclimatebootvpdcoef350km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wBclimatebootvpdcoef840km.csv',sep=""),
            sep=",",row.names=F)
write.table(Bclimatebootmagcoef350km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wBclimatebootmagcoef840km.csv',sep=""),
            sep=",",row.names=F)
write.table(Bclimatebootfreezecoef350km,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wBclimatebootfreezecoef840km.csv',sep=""),
            sep=",",row.names=F)
setwd(wname)
#
#
#
# Aggregate the results of the regression coefficients
# ----------------------------------------------------
slopemeanmindist <- data.frame(horizon=c("A","Ap","B"),
                               ppt=NA,
                               vpd=NA,
                               mag=NA,
                               freeze=NA)
interceptmeanmindist <- slopemeanmindist
slopesdmindist <- slopemeanmindist
interceptsdmindist <- slopemeanmindist
slopePvaluemindist <- slopemeanmindist
interceptPvaluemindist <- slopemeanmindist
Propspatstructrm <- slopemeanmindist # portion of the regressions with 
# spatial structure
bootAlist <- list(Aclimatebootpptcoef600km,
                  Aclimatebootvpdcoef600km,
                  Aclimatebootmagcoef600km,
                  Aclimatebootfreezecoef600km)
bootAplist <- list(Apclimatebootpptcoef220km,
                   Apclimatebootvpdcoef220km,
                   Apclimatebootmagcoef220km,
                   Apclimatebootfreezecoef220km)
bootBlist <- list(Bclimatebootpptcoef350km,
                  Bclimatebootvpdcoef350km,
                  Bclimatebootmagcoef350km,
                  Bclimatebootfreezecoef350km)
for (i in 1:length(bootAlist)) {
  slopemeanmindist[1,i+1] <- mean(bootAlist[[i]][,"slope"])
  interceptmeanmindist[1,i+1] <- mean(bootAlist[[i]][,"intercept"])
  slopesdmindist[1,i+1] <- sd(bootAlist[[i]][,"slope"])
  interceptsdmindist[1,i+1] <- sd(bootAlist[[i]][,"intercept"])
  slopePvaluemindist[1,i+1] <- t.test(bootAlist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindist[1,i+1] <- t.test(bootAlist[[i]][,"intercept"],
                                          mu=0)$p.value
  Propspatstructrm[1,i+1] <- length(which(bootAlist[[i]][,"Ipvalue"] > 
                                            0.05))/N
  slopemeanmindist[2,i+1] <- mean(bootAplist[[i]][,"slope"])
  interceptmeanmindist[2,i+1] <- mean(bootAplist[[i]][,"intercept"])
  slopesdmindist[2,i+1] <- sd(bootAplist[[i]][,"slope"])
  interceptsdmindist[2,i+1] <- sd(bootAplist[[i]][,"intercept"])
  slopePvaluemindist[2,i+1] <- t.test(bootAplist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindist[2,i+1] <- t.test(bootAplist[[i]][,"intercept"],
                                          mu=0)$p.value
  Propspatstructrm[2,i+1] <- length(which(bootAplist[[i]][,"Ipvalue"] > 
                                            0.05))/N
  slopemeanmindist[3,i+1] <- mean(bootBlist[[i]][,"slope"])
  interceptmeanmindist[3,i+1] <- mean(bootBlist[[i]][,"intercept"])
  slopesdmindist[3,i+1] <- sd(bootBlist[[i]][,"slope"])
  interceptsdmindist[3,i+1] <- sd(bootBlist[[i]][,"intercept"])
  slopePvaluemindist[3,i+1] <- t.test(bootBlist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindist[3,i+1] <- t.test(bootBlist[[i]][,"intercept"],
                                          mu=0)$p.value
  Propspatstructrm[3,i+1] <- length(which(bootBlist[[i]][,"Ipvalue"] > 
                                            0.05))/N
}
#
#
# Aggregate the results of the regression coefficients after accounting
# for the spatial structure columns
# ---------------------------------------------------------------------
slopemeanmindistspatstructrm <- data.frame(horizon=c("A","Ap","B"),
                                           ppt=NA,
                                           vpd=NA,
                                           mag=NA,
                                           freeze=NA)
interceptmeanmindistspatstructrm <- slopemeanmindist
slopesdmindistspatstructrm <- slopemeanmindist
interceptsdmindistspatstructrm <- slopemeanmindist
slopePvaluemindistspatstructrm <- slopemeanmindist
interceptPvaluemindistspatstructrm <- slopemeanmindist
bootAlist <- list(Aclimatebootpptcoef600km[which(
  Aclimatebootpptcoef600km$Ipvalue > 0.05)[1:1200],],
                  Aclimatebootvpdcoef600km[which(
                    Aclimatebootvpdcoef600km$Ipvalue > 0.05)[1:1200],],
                  Aclimatebootmagcoef600km[which(
                    Aclimatebootmagcoef600km$Ipvalue > 0.05)[1:1200],],
                  Aclimatebootfreezecoef600km[which(
                    Aclimatebootfreezecoef600km$Ipvalue > 0.05)[1:1200],])
bootAplist <- list(Apclimatebootpptcoef220km[which(
  Apclimatebootpptcoef220km$Ipvalue > 0.05)[1:1200],],
                   Apclimatebootvpdcoef220km[which(
                     Apclimatebootvpdcoef220km$Ipvalue > 0.05)[1:1200],],
                   Apclimatebootmagcoef220km[which(
                     Apclimatebootmagcoef220km$Ipvalue > 0.05)[1:1200],],
                   Apclimatebootfreezecoef220km[which(
                     Apclimatebootfreezecoef220km$Ipvalue > 0.05)[1:1200],])
bootBlist <- list(Bclimatebootpptcoef350km[which(
  Bclimatebootpptcoef350km$Ipvalue > 0.05)[1:1200],],
                  Bclimatebootvpdcoef350km[which(
                    Bclimatebootvpdcoef350km$Ipvalue > 0.05)[1:1200],],
                  Bclimatebootmagcoef350km[which(
                    Bclimatebootmagcoef350km$Ipvalue > 0.05)[1:1200],],
                  Bclimatebootfreezecoef350km[which(
                    Bclimatebootfreezecoef350km$Ipvalue > 0.05)[1:1200],])
for (i in 1:length(bootAlist)) {
  slopemeanmindistspatstructrm[1,i+1] <- mean(
    bootAlist[[i]][,"slope"])
  interceptmeanmindistspatstructrm[1,i+1] <- mean(
    bootAlist[[i]][,"intercept"])
  slopesdmindistspatstructrm[1,i+1] <- sd(
    bootAlist[[i]][,"slope"])
  interceptsdmindistspatstructrm[1,i+1] <- sd(
    bootAlist[[i]][,"intercept"])
  slopePvaluemindistspatstructrm[1,i+1] <- t.test(
    bootAlist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindistspatstructrm[1,i+1] <- t.test(
    bootAlist[[i]][,"intercept"],mu=0)$p.value
  #
  slopemeanmindistspatstructrm[2,i+1] <- mean(
    bootAplist[[i]][,"slope"])
  interceptmeanmindistspatstructrm[2,i+1] <- mean(
    bootAplist[[i]][,"intercept"])
  slopesdmindistspatstructrm[2,i+1] <- sd(
    bootAplist[[i]][,"slope"])
  interceptsdmindistspatstructrm[2,i+1] <- sd(
    bootAplist[[i]][,"intercept"])
  slopePvaluemindistspatstructrm[2,i+1] <- t.test(
    bootAplist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindistspatstructrm[2,i+1] <- t.test(
    bootAplist[[i]][,"intercept"],mu=0)$p.value
  #
  slopemeanmindistspatstructrm[3,i+1] <- mean(
    bootBlist[[i]][,"slope"])
  interceptmeanmindistspatstructrm[3,i+1] <- mean(
    bootBlist[[i]][,"intercept"])
  slopesdmindistspatstructrm[3,i+1] <- sd(
    bootBlist[[i]][,"slope"])
  interceptsdmindistspatstructrm[3,i+1] <- sd(
    bootBlist[[i]][,"intercept"])
  slopePvaluemindistspatstructrm[3,i+1] <- t.test(
    bootBlist[[i]][,"slope"],mu=0)$p.value
  interceptPvaluemindistspatstructrm[3,i+1] <- t.test(
    bootBlist[[i]][,"intercept"],mu=0)$p.value
}
#
#
# Note: If this works, need to make plots of using these values and plot
# them to a new Figure 2 like figure.
# [insert code here]
#
#
# Output the summary of the results of the bootstrapped aggregated data 
# separated by the appropriate minimum distance
# ---------------------------------------------------------------------
setwd(oname)
# w = weighted regression
write.table(slopemeanmindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wslopemeanmindist.csv',sep=""),
            sep=",",row.names=F)
write.table(interceptmeanmindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_winterceptmeanmindist.csv',sep=""),
            sep=",",row.names=F)
write.table(slopesdmindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wslopesdmindist.csv',sep=""),
            sep=",",row.names=F)
write.table(interceptsdmindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_winterceptsdmindist.csv',sep=""),
            sep=",",row.names=F)
write.table(slopePvaluemindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wslopePvaluemindist.csv',sep=""),
            sep=",",row.names=F)
write.table(interceptPvaluemindist,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_winterceptPvaluemindist.csv',sep=""),
            sep=",",row.names=F)
write.table(Propspatstructrm,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wPropspatstructrm.csv',sep=""),
            sep=",",row.names=F)
#
#
write.table(slopemeanmindistspatstructrm,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wslopemeanmindistspatstructrm.csv',sep=""),
            sep=",",row.names=F)
write.table(interceptmeanmindistspatstructrm,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_winterceptmeanmindistspatstructrm.csv',sep=""),
            sep=",",row.names=F)
write.table(slopesdmindistspatstructrm,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wslopesdmindistspatstructrm.csv',sep=""),
            sep=",",row.names=F)
write.table(interceptsdmindistspatstructrm,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_winterceptsdmindistspatstructrm.csv',sep=""),
            sep=",",row.names=F)
write.table(slopePvaluemindistspatstructrm,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_wslopePvaluemindistspatstructrm.csv',sep=""),
            sep=",",row.names=F)
write.table(interceptPvaluemindistspatstructrm,
            file=paste(gsub("-", "", 
                            as.character(Sys.Date())),
                       '_winterceptPvaluemindistspatstructrm.csv',sep=""),
            sep=",",row.names=F)
setwd(wname)
#
# 
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Recreate Figure 2 from the revised submission (modified from the "R log -
# new calculations.R" script)
# ---------------------------------------------------------------------------
setwd(fname)
# Open the plotting region
figuretype <- "eps" # either pdf or eps
if (figuretype == "pdf") {
  pdf(file=paste(gsub("-", "", 
                      as.character(Sys.Date())),
                 '_ClimateStationResidualEP.pdf',sep=""),
      height=7,width=8.5)
} else { # .eps
  postscript(file=paste(gsub("-", "", 
                             as.character(Sys.Date())),
                        '_ClimateStationResidualEP.eps',sep=""),
             height=7,width=8.5)
}
# Set graphical parameters for a multipanel plot
par(mfrow=c(3,4),omi=c(1,0.5,1,0.5),mar=c(3.2,3.2,1,1),mgp=c(2.2,1,0))
#ylimhold <- c(-0.22,.18) # defines the ylim scale for all the plots below
ylimhold <- c(-0.18,.18) # defines the ylim scale for all the plots below
#
# Set the xlim for each plot (determined after running the outlier 
# detection)
pptxlim <- c(50,1750)
vpdxlim <- c(4,46)
magxlim <- c(2,18)
freezexlim <- c(0,290)
#
# Save the summary coefficients from all the work above (resampling 
# methods)
slopesummary <- slopemeanmindistspatstructrm # slopes
intsummary <- interceptmeanmindistspatstructrm # intercepts
slopesdsummary <- slopesdmindistspatstructrm # slope standard deviations
# based on the regression slopes from the 1200 resampled subsets
slopesesummary <- slopesdsummary # standard errors
slopesesummary$ppt <- slopesdsummary$ppt/sqrt(1200)
slopesesummary$vpd <- slopesdsummary$vpd/sqrt(1200)
slopesesummary$mag <- slopesdsummary$mag/sqrt(1200)
slopesesummary$freeze <- slopesdsummary$freeze/sqrt(1200)
#
# Critical t value
alpha <- 0.05 # alpha-level for a 95% confidence interval (alpha is used
# far above in the calculation of the Poisson distribution summary from
# the weather station data)
tc <- qt(1-alpha/2,1200-1) # 95% (1-0.05/2)
# 
# A horizons
# ----------
colorhold <- "gainsboro"
hor <- "A"
#
# Residual effective porosity vs MAP
dhold <- Aclimatex2
summarycol <- "ppt"
pred <- "ppt"
dhold <- dhold[which(dhold$resepxpptout == FALSE),]
nApptx <- nrow(dhold)
plot(resepx~ppt,dhold,pch=21,bg=colorhold,
     xlab="Mean annual precipitation (mm)",
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=pptxlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# Residual effective porosity vs vapor pressure deficit
dhold <- Aclimatex2
summarycol <- "vpd"
pred <- "prismjjavpd"
dhold <- dhold[which(dhold$resepxprismjjavpdout == FALSE),]
nAvpdx <- nrow(dhold)
plot(resepx~prismjjavpd,dhold,pch=21,bg=colorhold,
     xlab="Mean max JJA vapor pressure deficit (hPa)",
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=vpdxlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# Residual effective porosity vs precipitation magnitude
dhold <- Aclimatex2
summarycol <- "mag"
pred <- "magnitude"
dhold <- dhold[which(dhold$resepxmagnitudeout == FALSE),]
nAmagx <- nrow(dhold)
plot(resepx~magnitude,dhold,pch=21,bg=colorhold,
     xlab="Precipitation event magnitude (mm)",
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=magxlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# Residual effective porosity vs freeze frequency
dhold <- Aclimatex2
summarycol <- "freeze"
pred <- "mintempfreezedays"
dhold <- dhold[which(dhold$resepxmintempfreezedaysout == FALSE),]
nAfreezex <- nrow(dhold)
plot(resepx~mintempfreezedays,dhold,pch=21,bg=colorhold,
     xlab=expression(paste("Freezing frequency (N yr"^{-1},")",
                           sep="")),
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=freezexlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# Ap horizons
# -----------
colorhold <- "darkseagreen3"
hor <- "Ap"
#
# Residual effective porosity vs MAP
dhold <- Apclimatex2
summarycol <- "ppt"
pred <- "ppt"
dhold <- dhold[which(dhold$resepxpptout == FALSE),]
nAppptx <- nrow(dhold)
plot(resepx~ppt,dhold,pch=21,bg=colorhold,
     xlab="Mean annual precipitation (mm)",
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=pptxlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# Residual effective porosity vs vapor pressure deficit
dhold <- Apclimatex2
summarycol <- "vpd"
pred <- "prismjjavpd"
dhold <- dhold[which(dhold$resepxprismjjavpdout == FALSE),]
nApvpdx <- nrow(dhold)
plot(resepx~prismjjavpd,dhold,pch=21,bg=colorhold,
     xlab="Mean max JJA vapor pressure deficit (hPa)",
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=vpdxlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# Residual effective porosity vs precipitation magnitude
dhold <- Apclimatex2
summarycol <- "mag"
pred <- "magnitude"
dhold <- dhold[which(dhold$resepxmagnitudeout == FALSE),]
nApmagx <- nrow(dhold)
plot(resepx~magnitude,dhold,pch=21,bg=colorhold,
     xlab="Precipitation event magnitude (mm)",
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=magxlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# Residual effective porosity vs freeze frequency
dhold <- Apclimatex2
summarycol <- "freeze"
pred <- "mintempfreezedays"
dhold <- dhold[which(dhold$resepxmintempfreezedaysout == FALSE),]
nApfreezex <- nrow(dhold)
plot(resepx~mintempfreezedays,dhold,pch=21,bg=colorhold,
     xlab=expression(paste("Freezing frequency (N yr"^{-1},")",
                           sep="")),
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=freezexlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# B horizons
# ----------
colorhold <- "lightcoral"
dhold <- Bclimatex2
hor <- "B"
#
# Residual effective porosity vs MAP
dhold <- Bclimatex2
summarycol <- "ppt"
pred <- "ppt"
dhold <- dhold[which(dhold$resepxpptout == FALSE),]
nBpptx <- nrow(dhold)
plot(resepx~ppt,dhold,pch=21,bg=colorhold,
     xlab="Mean annual precipitation (mm)",
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=pptxlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# Residual effective porosity vs vapor pressure deficit
dhold <- Bclimatex2
summarycol <- "vpd"
pred <- "prismjjavpd"
dhold <- dhold[which(dhold$resepxprismjjavpdout == FALSE),]
nBvpdx <- nrow(dhold)
plot(resepx~prismjjavpd,dhold,pch=21,bg=colorhold,
     xlab="Mean max JJA vapor pressure deficit (hPa)",
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=vpdxlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# Residual effective porosity vs precipitation magnitude
dhold <- Bclimatex2
summarycol <- "mag"
pred <- "magnitude"
dhold <- dhold[which(dhold$resepxmagnitudeout == FALSE),]
nBmagx <- nrow(dhold)
plot(resepx~magnitude,dhold,pch=21,bg=colorhold,
     xlab="Precipitation event magnitude (mm)",
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=magxlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# Residual effective porosity vs freeze frequency
dhold <- Bclimatex2
summarycol <- "freeze"
pred <- "mintempfreezedays"
dhold <- dhold[which(dhold$resepxmintempfreezedaysout == FALSE),]
nBfreezex <- nrow(dhold)
plot(resepx~mintempfreezedays,dhold,pch=21,bg=colorhold,
     xlab=expression(paste("Freezing frequency (N yr"^{-1},")",
                           sep="")),
     ylab="Residual effective porosity",
     ylim=ylimhold,
     xlim=freezexlim,
     cex=1.2)
m <- slopesummary[which(slopesummary$horizon==hor),summarycol]
mse <- slopesesummary[which(slopesummary$horizon==hor),summarycol]
b <- intsummary[which(slopesummary$horizon==hor),summarycol]
x <- seq(-2000,2000,length.out=10000)
y <- m*x + b
lines(y~x,col="black",lwd=1.5,lty=1) # bootstrapped OLS
lCI <- m - tc*mse # lower confidence interval
uCI <- m + tc*mse # upper confidence interval
centr <- mean(dhold[,pred]) # centroid (x)
centr <- c(centr, m*centr + b) # centroid (x and y) that falls on the
# line
# Lower slope defined by the 95% CI:
deltay <- (lCI*centr[1]+b) - centr[2]
y <- lCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
# Upper slope defined by the 95% CI:
deltay <- (uCI*centr[1]+b) - centr[2]
y <- uCI*x+b
lines(I(y-deltay)~x,col="blue",lwd=1.5,lty=2)
lines(c(-10000,10000),c(0,0),lty=2)
#
# Close the device
dev.off()
setwd(wname)
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Recreate Figure 3 from the revised submission
# ---------------------------------------------------------------------------
#
# Calculates the change in precipitation and effective porosity for the 
# Northern Great Plains, Southern Great Plains, Basin and Range,
# Central Lowlands, Southeast Coastal Plain, and Pacific Northwest
# ---------------------------------------------------------------------
#
# Calculates the indices of the climate dataframe that correspond to
# the relevant regions (NOTE!!! I use climate instead of climatex
# since I use the original EP values below and climatex has lost
# them somewhere on the way)
ngp <- which(((climate$lat < 48.63) & (climate$lat > 43.32)) &
               ((climate$long < -92.44) & (climate$long > -105.36)))
sgp <- which(((climate$lat < 35.6) & (climate$lat > 30.06)) &
               ((climate$long < -95.16) & (climate$long > -104.66)))
br <- which(((climate$lat < 41.07) & (climate$lat > 36.2)) &
              ((climate$long < -110.55) & (climate$long > -119.51)))
clow <- which(((climate$lat < 41.27) & (climate$lat > 37.33)) &
                ((climate$long < -83.39) & (climate$long > -90.51)))
scp <- which(((climate$lat < 34.41) & (climate$lat > 30.86)) &
               ((climate$long < -84.09) & (climate$long > -92.97)))
pn <- which(((climate$lat < 47.9) & (climate$lat > 44.74)) &
              ((climate$long < -115.03) & (climate$long > -120.92)))
#
# Creates new dataframes of A horizon properties near climate stations
# within the various regions
cAngp <- climate[ngp,]
cAsgp <- climate[sgp,]
cAbr <- climate[br,]
cAclow <- climate[clow,]
cAscp <- climate[scp,]
cApn <- climate[pn,]
#
# Predicts the Ks using the current EP
Kmethod <- "Rawls" # option to predict Ksat using either Rawls et al. 
# (1998) or Han et al. (1998), "Han"
if (Kmethod == "Rawls") {
  hold <- cAngp
  cAngp$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],
                       unit="cmd",method=Kmethod) # cm/d...
  hold <- cAsgp
  cAsgp$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],
                       unit="cmd",method=Kmethod)
  hold <- cAbr
  cAbr$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],
                      unit="cmd",method=Kmethod)
  hold <- cAclow
  cAclow$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],
                        unit="cmd",method=Kmethod)
  hold <- cAscp
  cAscp$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],
                       unit="cmd",method=Kmethod)
  hold <- cApn
  cApn$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],
                      unit="cmd",method=Kmethod)
}
if (Kmethod == "Han") {
  hold <- cAngp
  cAngp$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],hold[,"Atp"],
                       unit="cmd",method=Kmethod) # cm/d...
  hold <- cAsgp
  cAsgp$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],hold[,"Atp"],
                       unit="cmd",method=Kmethod)
  hold <- cAbr
  cAbr$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],hold[,"Atp"],
                      unit="cmd",method=Kmethod)
  hold <- cAclow
  cAclow$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],hold[,"Atp"],
                        unit="cmd",method=Kmethod)
  hold <- cAscp
  cAscp$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],hold[,"Atp"],
                       unit="cmd",method=Kmethod)
  hold <- cApn
  cApn$Ks <- Ksatpred(hold[,"Aep"],hold[,"Afc"],hold[,"Awp"],hold[,"Atp"],
                      unit="cmd",method=Kmethod)
}
#
# Obtains the slope from the resampled OLS regressions for the A horizons
Aslope <- slopemeanmindistspatstructrm[which(
  slopemeanmindistspatstructrm=="A"),"ppt"]
#
# Predict the change in EP as a result of the change in precipitation
# at the end of the century; this is the same as predicting the change
# in residual EP; also predict the new EP
hold <- cAngp
cAngp$deltaAep <- Aslope*(hold$cesmppt-hold$prismppt)
cAngp$predAep <- cAngp$Aep+cAngp$deltaAep
hold <- cAsgp
cAsgp$deltaAep <- Aslope*(hold$cesmppt-hold$prismppt)
cAsgp$predAep <- cAsgp$Aep+cAsgp$deltaAep
hold <- cAbr
cAbr$deltaAep <- Aslope*(hold$cesmppt-hold$prismppt)
cAbr$predAep <- cAbr$Aep+cAbr$deltaAep
hold <- cAclow
cAclow$deltaAep <- Aslope*(hold$cesmppt-hold$prismppt)
cAclow$predAep <- cAclow$Aep+cAclow$deltaAep
hold <- cAscp
cAscp$deltaAep <- Aslope*(hold$cesmppt-hold$prismppt)
cAscp$predAep <- cAscp$Aep+cAscp$deltaAep
hold <- cApn
cApn$deltaAep <- Aslope*(hold$cesmppt-hold$prismppt)
cApn$predAep <- cApn$Aep+cApn$deltaAep
#
# Predict the Ks that results from the future EP which changes as
# a result of precipitation changes
if (Kmethod == "Rawls") {
  hold <- cAngp
  cAngp$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],hold[,"Awp"],
                           unit="cmd",method=Kmethod) # cm/d...
  hold <- cAsgp
  cAsgp$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],hold[,"Awp"],
                           unit="cmd",method=Kmethod)
  hold <- cAbr
  cAbr$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],hold[,"Awp"],
                          unit="cmd",method=Kmethod)
  hold <- cAclow
  cAclow$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],hold[,"Awp"],
                            unit="cmd",method=Kmethod)
  hold <- cAscp
  cAscp$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],hold[,"Awp"],
                           unit="cmd",method=Kmethod)
  hold <- cApn
  cApn$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],hold[,"Awp"],
                          unit="cmd",method=Kmethod)
}
if (Kmethod == "Han") {
  hold <- cAngp
  cAngp$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],
                           hold[,"Awp"],hold[,"Atp"],
                           unit="cmd",method=Kmethod) # cm/d...
  hold <- cAsgp
  cAsgp$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],
                           hold[,"Awp"],hold[,"Atp"],
                           unit="cmd",method=Kmethod)
  hold <- cAbr
  cAbr$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],
                          hold[,"Awp"],hold[,"Atp"],
                          unit="cmd",method=Kmethod)
  hold <- cAclow
  cAclow$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],
                            hold[,"Awp"],hold[,"Atp"],
                            unit="cmd",method=Kmethod)
  hold <- cAscp
  cAscp$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],
                           hold[,"Awp"],hold[,"Atp"],
                           unit="cmd",method=Kmethod)
  hold <- cApn
  cApn$predKs <- Ksatpred(hold[,"predAep"],hold[,"Afc"],
                          hold[,"Awp"],hold[,"Atp"],
                          unit="cmd",method=Kmethod)
}
#
# Calculate the change in conductivity and the percent change
hold <- cAngp
cAngp$deltaKs <- hold$predKs-hold$Ks # change
cAngp$deltaKsper <- cAngp$deltaKs/cAngp$Ks*100 # percent change
cAngp$region <- "Northern Great Plains"
hold <- cAsgp
cAsgp$deltaKs <- hold$predKs-hold$Ks # change
cAsgp$deltaKsper <- cAsgp$deltaKs/cAsgp$Ks*100 # percent change
cAsgp$region <- "Southern Great Plains"
hold <- cAbr
cAbr$deltaKs <- hold$predKs-hold$Ks # change
cAbr$deltaKsper <- cAbr$deltaKs/cAbr$Ks*100 # percent change
cAbr$region <- "Basin and Range"
hold <- cAclow
cAclow$deltaKs <- hold$predKs-hold$Ks # change
cAclow$deltaKsper <- cAclow$deltaKs/cAclow$Ks*100 # percent change
cAclow$region <- "Central Lowlands"
hold <- cAscp
cAscp$deltaKs <- hold$predKs-hold$Ks # change
cAscp$deltaKsper <- cAscp$deltaKs/cAscp$Ks*100 # percent change
cAscp$region <- "Southeast Coastal Plain"
hold <- cApn
cApn$deltaKs <- hold$predKs-hold$Ks # change
cApn$deltaKsper <- cApn$deltaKs/cApn$Ks*100 # percent change
cApn$region <- "Pacific Northwest"
#
#
# Combine the individual climate data.frames together with the 
# appropriate name
# ------------------------------------------------------------
regdata <- rbind(cAngp,cAsgp,cAbr,cAclow,cAscp,cApn)
#
#
# Plot of histograms showing predicted changes in Ks for the various
# regions
# ------------------------------------------------------------------
histhold <- regdata[which(regdata$region=="Northern Great Plains"),]
hngp <- hist(histhold$deltaKsper,breaks=6)
histhold <- regdata[which(regdata$region=="Southern Great Plains"),]
hsgp <- hist(histhold$deltaKsper,breaks=10)
histhold <- regdata[which(regdata$region=="Basin and Range"),]
hbr <- hist(histhold$deltaKsper,breaks=10)
histhold <- regdata[which(regdata$region=="Central Lowlands"),]
hclow <- hist(histhold$deltaKsper,breaks=6)
histhold <- regdata[which(regdata$region=="Southeast Coastal Plain"),]
hscp <- hist(histhold$deltaKsper,breaks=10)
histhold <- regdata[which(regdata$region=="Pacific Northwest"),]
hpn <- hist(histhold$deltaKsper,breaks=10)
#
# Open the plotting region
setwd(fname)
figuretype <- "eps" # either pdf or eps
if (figuretype == "pdf") {
  pdf(file=paste(gsub("-", "", 
                      as.character(Sys.Date())),
                 '_Ksxhistograms.pdf',sep=""),
      height=7,width=7)
} else { # .eps
  postscript(file=paste(gsub("-", "", 
                             as.character(Sys.Date())),
                        '_Ksxhistograms.eps',sep=""),
             height=7,width=7)
}
xlimhold <- c(-60,60)
ylimhold <- c(0,0.06)
# Set graphical parameters for a multipanel plot
par(mfrow=c(3,3),omi=c(1,0.5,1,0.5),mar=c(3.2,3.2,1,1),mgp=c(2.2,1,0))
negcolor <- "red2"
poscolor <- "royalblue4"
#
cuts <- cut(hngp$breaks, c(-Inf,-0.5,Inf))
cuts <- as.character(cuts)
cuts[which(cuts=="(-Inf,-0.5]")] <- negcolor
cuts[which(cuts=="(-0.5, Inf]")] <- poscolor
plot(hngp, col=cuts,freq=FALSE,
     xlab=expression(paste(Delta,"K"["s"]," (%)",sep="")),
     ylab="Density",
     main="Northern Great Plains",
     xlim=xlimhold,
     ylim=ylimhold,
     cex=1.2)
box()
cuts <- cut(hsgp$breaks, c(-Inf,-0.5,Inf))
cuts <- as.character(cuts)
cuts[which(cuts=="(-Inf,-0.5]")] <- negcolor
cuts[which(cuts=="(-0.5, Inf]")] <- poscolor
plot(hsgp, col=cuts,freq=FALSE,
     xlab=expression(paste(Delta,"K"["s"]," (%)",sep="")),
     ylab="Density",
     main="Southern Great Plains",
     xlim=xlimhold,
     ylim=ylimhold,
     cex=1.2)
box()
cuts <- cut(hbr$breaks, c(-Inf,-0.5,Inf))
cuts <- as.character(cuts)
cuts[which(cuts=="(-Inf,-0.5]")] <- negcolor
cuts[which(cuts=="(-0.5, Inf]")] <- poscolor
plot(hbr, col=cuts,freq=FALSE,
     xlab=expression(paste(Delta,"K"["s"]," (%)",sep="")),
     ylab="Density",
     main="Basin and Range",
     xlim=xlimhold,
     ylim=ylimhold,
     cex=1.2)
box()
cuts <- cut(hclow$breaks, c(-Inf,-0.5,Inf))
cuts <- as.character(cuts)
cuts[which(cuts=="(-Inf,-0.5]")] <- negcolor
cuts[which(cuts=="(-0.5, Inf]")] <- poscolor
plot(hclow, col=cuts,freq=FALSE,
     xlab=expression(paste(Delta,"K"["s"]," (%)",sep="")),
     ylab="Density",
     main="Central Lowlands",
     xlim=xlimhold,
     ylim=ylimhold,
     cex=1.2)
box()
cuts <- cut(hscp$breaks, c(-Inf,-0.5,Inf))
cuts <- as.character(cuts)
cuts[which(cuts=="(-Inf,-0.5]")] <- negcolor
cuts[which(cuts=="(-0.5, Inf]")] <- poscolor
plot(hscp, col=cuts,freq=FALSE,
     xlab=expression(paste(Delta,"K"["s"]," (%)",sep="")),
     ylab="Density",
     main="Southeast Coastal Plain",
     xlim=xlimhold,
     ylim=ylimhold,
     cex=1.2)
box()
cuts <- cut(hpn$breaks, c(-Inf,-0.5,Inf))
cuts <- as.character(cuts)
cuts[which(cuts=="(-Inf,-0.5]")] <- negcolor
cuts[which(cuts=="(-0.5, Inf]")] <- poscolor
plot(hpn, col=cuts,freq=FALSE,
     xlab=expression(paste(Delta,"K"["s"]," (%)",sep="")),
     ylab="Density",
     main="Pacific Northwest",
     xlim=xlimhold,
     ylim=ylimhold,
     cex=1.2)
box()
dev.off()
setwd(wname)
#
# Some data for the plots:
aggregate(deltaKsper ~ region, regdata, FUN = each(mean,min,max))
# ---------------------------------------------------------------------------