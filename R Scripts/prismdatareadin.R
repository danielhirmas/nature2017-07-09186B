# ---------------------------------------------------------------------------
# Title: R log for to read in the PRISM Data
# Author: D.R. Hirmas
# Date: 2/1/2018
# ---------------------------------------------------------------------------
#
#
# ---------------------------------------------------------------------------
# Read in the data
# ---------------------------------------------------------------------------
# Change the working directories and read in the appropriate files
dname <- paste(dname,"PRISM Climate Data for the US",sep="/")
#
# Mean Annual Precipitation
# -------------------------
setwd(paste(dname,"PRISM_ppt_30yr_normal_4kmM2_all_asc",sep="/"))
prismppt <- read.table("PRISM_ppt_30yr_normal_4kmM2_annual_asc.asc",
                  skip=6,
                  sep=" ",
                  header=F)
prismppt <- prismppt[,-1] # removes the space at the beginning of each row in the data
prismppt <- as.matrix(prismppt) # saves as a matrix
colnames(prismppt) <- NULL # removes the column names
prismppt[which(prismppt==-9999)] <- NA # removes the NA values
#
# Mean Annual Temperature
# -----------------------
setwd(paste(dname,"PRISM_tmean_30yr_normal_4kmM2_all_asc",sep="/"))
prismtemp <- read.table("PRISM_tmean_30yr_normal_4kmM2_annual_asc.asc",
                  skip=6,
                  sep=" ",
                  header=F)
prismtemp <- prismtemp[,-1] # removes the space at the beginning of each row in the data
prismtemp <- as.matrix(prismtemp) # saves as a matrix
colnames(prismtemp) <- NULL # removes the column names
prismtemp[which(prismtemp==-9999)] <- NA # removes the NA values
#
# Mean Dewpoint Temperature
# -------------------------
setwd(paste(dname,"PRISM_tdmean_30yr_normal_4kmM2_all_asc",sep="/"))
prismdewtemp <- read.table("PRISM_tdmean_30yr_normal_4kmM2_annual_asc.asc",
                   skip=6,
                   sep=" ",
                   header=F)
prismdewtemp <- prismdewtemp[,-1] # removes the space at the beginning of each row in the data
prismdewtemp <- as.matrix(prismdewtemp) # saves as a matrix
colnames(prismdewtemp) <- NULL # removes the column names
prismdewtemp[which(prismdewtemp==-9999)] <- NA # removes the NA values
#
# Mean Annual Minumum Temperature
# -------------------------------
setwd(paste(dname,"PRISM_tmin_30yr_normal_4kmM2_all_asc",sep="/"))
prismmintemp <- read.table("PRISM_tmin_30yr_normal_4kmM2_annual_asc.asc",
                      skip=6,
                      sep=" ",
                      header=F)
prismmintemp <- prismmintemp[,-1] # removes the space at the beginning of each row in the data
prismmintemp <- as.matrix(prismmintemp) # saves as a matrix
colnames(prismmintemp) <- NULL # removes the column names
prismmintemp[which(prismmintemp==-9999)] <- NA # removes the NA values
#
# Mean Annual Maximum Temperature
# -------------------------------
setwd(paste(dname,"PRISM_tmax_30yr_normal_4kmM2_all_asc",sep="/"))
prismmaxtemp <- read.table("PRISM_tmax_30yr_normal_4kmM2_annual_asc.asc",
                      skip=6,
                      sep=" ",
                      header=F)
prismmaxtemp <- prismmaxtemp[,-1] # removes the space at the beginning of each row in the data
prismmaxtemp <- as.matrix(prismmaxtemp) # saves as a matrix
colnames(prismmaxtemp) <- NULL # removes the column names
prismmaxtemp[which(prismmaxtemp==-9999)] <- NA # removes the NA values
#
# Mean Minumum Vapor Pressure Deficit
# -----------------------------------
setwd(paste(dname,"PRISM_vpdmin_30yr_normal_4kmM2_all_asc",sep="/"))
prismminvpd <- read.table("PRISM_vpdmin_30yr_normal_4kmM2_annual_asc.asc",
                      skip=6,
                      sep=" ",
                      header=F)
prismminvpd <- prismminvpd[,-1] # removes the space at the beginning of each row in the data
prismminvpd <- as.matrix(prismminvpd) # saves as a matrix
colnames(prismminvpd) <- NULL # removes the column names
prismminvpd[which(prismminvpd==-9999)] <- NA # removes the NA values
#
# Mean Maximum Vapor Pressure Deficit
# -----------------------------------
setwd(paste(dname,"PRISM_vpdmax_30yr_normal_4kmM2_all_asc",sep="/"))
prismmaxvpd <- read.table("PRISM_vpdmax_30yr_normal_4kmM2_annual_asc.asc",
                     skip=6,
                     sep=" ",
                     header=F)
prismmaxvpd <- prismmaxvpd[,-1] # removes the space at the beginning of each row in the data
prismmaxvpd <- as.matrix(prismmaxvpd) # saves as a matrix
colnames(prismmaxvpd) <- NULL # removes the column names
prismmaxvpd[which(prismmaxvpd==-9999)] <- NA # removes the NA values
#
# Mean Maximum June Vapor Pressure Deficit
# ----------------------------------------
setwd(paste(dname,"PRISM_vpdmax_30yr_normal_4kmM2_06_asc",sep="/"))
prismjunemaxvpd <- read.table("PRISM_vpdmax_30yr_normal_4kmM2_06_asc.asc",
                          skip=6,
                          sep=" ",
                          header=F)
prismjunemaxvpd <- prismjunemaxvpd[,-1] # removes the space at the beginning of each row in the data
prismjunemaxvpd <- as.matrix(prismjunemaxvpd) # saves as a matrix
colnames(prismjunemaxvpd) <- NULL # removes the column names
prismjunemaxvpd[which(prismjunemaxvpd==-9999)] <- NA # removes the NA values
#
# Mean Maximum July Vapor Pressure Deficit
# ----------------------------------------
setwd(paste(dname,"PRISM_vpdmax_30yr_normal_4kmM2_07_asc",sep="/"))
prismjulymaxvpd <- read.table("PRISM_vpdmax_30yr_normal_4kmM2_07_asc.asc",
                              skip=6,
                              sep=" ",
                              header=F)
prismjulymaxvpd <- prismjulymaxvpd[,-1] # removes the space at the beginning of each row in the data
prismjulymaxvpd <- as.matrix(prismjulymaxvpd) # saves as a matrix
colnames(prismjulymaxvpd) <- NULL # removes the column names
prismjulymaxvpd[which(prismjulymaxvpd==-9999)] <- NA # removes the NA values
#
# Mean Maximum August Vapor Pressure Deficit
# ------------------------------------------
setwd(paste(dname,"PRISM_vpdmax_30yr_normal_4kmM2_08_asc",sep="/"))
prismaugustmaxvpd <- read.table("PRISM_vpdmax_30yr_normal_4kmM2_08_asc.asc",
                              skip=6,
                              sep=" ",
                              header=F)
prismaugustmaxvpd <- prismaugustmaxvpd[,-1] # removes the space at the beginning of each row in the data
prismaugustmaxvpd <- as.matrix(prismaugustmaxvpd) # saves as a matrix
colnames(prismaugustmaxvpd) <- NULL # removes the column names
prismaugustmaxvpd[which(prismaugustmaxvpd==-9999)] <- NA # removes the NA values
#
# Mean Maximum Daily June July and August Vapor Pressure Deficit
# --------------------------------------------------------------
prismjjavpd <- (prismjunemaxvpd+prismjulymaxvpd+prismaugustmaxvpd)/3
#
# Elevation
# ---------
setwd(paste(dname,"PRISM_us_dem_4km_asc",sep="/"))
elev <- read.table("PRISM_us_dem_4km_asc.asc",
                     skip=6,
                     sep=" ",
                     header=F)
elev <- elev[,-1] # removes the space at the beginning of each row in the data
elev <- as.matrix(elev) # saves as a matrix
colnames(elev) <- NULL # removes the column names
elev[which(elev==-9999)] <- NA # removes the NA values
#
# Metadata
# --------
# This information is obtained in the first 6 rows of each PRISM data file
pncols <- 1405
pnrows <- 621
xllcorner <- -125.020833333333
yllcorner <- 24.062499999795
cellsize <- 0.041666666667
colref <- seq(xllcorner,xllcorner+pncols*cellsize,length.out=pncols)
# refers to longitude of each column
rowref <- rev(seq(yllcorner,yllcorner+pnrows*cellsize,length.out=pnrows))
# refers to lattitude of each row

# Clean up
# --------
rm(dname)
# ---------------------------------------------------------------------------

