# ---------------------------------------------------------------------------
# Title: R log for to read in the CESM Data
# Author: D.R. Hirmas
# Date: 2/1/2018
# ---------------------------------------------------------------------------
#
#
# ---------------------------------------------------------------------------
# Read in the data
# ---------------------------------------------------------------------------
# Need one of the functions defined in functions R file below
setwd(funname)
source("20180612_EPProject_DataAnalysis_Functions.R")
#
# Change the working directories and read in the appropriate files
dname <- paste(dname,"CESM Data",sep="/")
#
# Mean Annual Precipitation for the RCP 6 scenario
# ------------------------------------------------
setwd(dname)
cesm <- load('CESM_2081_2100_rcp6.Rdata') # in mm
#
# Average years 2081-2100 for output
# ----------------------------------
sumppt <- annualppt[,,1] # 2081
for (i in 2:20) sumppt <- sumppt + annualppt[,,i] # goes from 2082-2100
cesmppt <- sumppt/20 # 2081-2100 mean
cesmppt <- counterclockwise(cesmppt) # properly setsup the matrix for the 
# lat and long variables below
#
# Metadata
# --------
# annualppt: size = 288 192 20
# annualtemp: size = 288 192 20
# lat: size = 192
# lon: size = 288
# These are global annual total ppt and average temp.
# The 20 are the years 1 = 2081, 20 = 2100.
lon[which(lon>180)] <- lon[which(lon>180)]-360 # fixes the longitude
# of the CESM data
ccolref <- lon
# refers to longitude of each column
crowref <- rev(lat)
# refers to lattitude of each row
# the "c" refers to CESM to distinguish it from the PRISM colref and rowref

# Clean up
# --------
rm(dname)
# ---------------------------------------------------------------------------

