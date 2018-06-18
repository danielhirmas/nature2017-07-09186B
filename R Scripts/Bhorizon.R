# ---------------------------------------------------------------------------
# Examine the relationship between climate and ep for B horizons (see
# the selection of B horizons below) that have their mid horizon depth 
# between 25 and 100 cm
# ---------------------------------------------------------------------------
setwd(fname)
# Select the B horizons
mtlt <- "BHorizonswithMidpointDepthsbetween25-100cm"
mindepth <- 25 # selects horizons that are less than or equal to this depth 
# (cm)
maxdepth <- 100 # selects horizons that are greater than or equal to this 
# depth (cm)
horizonhold <- c("Bt1","Bw1","Bw","Bt","Bt1","B","B1","Bw2","B2","Btg1",
                 "Bg1","Bt11","Bg","B11","Bss1","Bt21","Bt3","B12","Bwg1",
                 "Bss","Bg2","Btss","Bssg1","Bw11","B13","Bw3","Bg11",
                 "Bw12","B23","Bgss1","Bc1","Btm","Bwg","Bwc","Bwc1",
                 "Btc","Bt13","Btg12","Bssg2")
require(fields)
hhold <- NULL
for (i in 1:nrow(dat)) {
  if ((length(grep("B",dat$hzn_master[i])) != 0) & (
    is.na(dat$hzn_top[i])==F)) {
    if ((length(which(horizonhold==dat$hzn_desgn[i])) != 0) &
          ((dat$mid[i]>=mindepth) & (dat$mid[i]<=maxdepth))) {
      hhold <- c(hhold,i)
    }
  }
}
noB <- length(hhold)
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
dB <- d # saves the data frame for this selection
#
# Save the lowest point in Florida for removal below (i.e., FLremove==T
# is set above)
FLmin <- min(d$lat[which(d$state_admindiv_code=="FL")])
#
#
# Interpolated maps
# _________________
# Porosities
# Set the scale after running both the A and the B horizons
minvalueep <- 0.09; maxvalueep <- 0.3
minvaluewhc <- 0.06; maxvaluewhc <- 0.22
pdf(file=paste(mtlt,"Interplotated Maps.pdf",sep="_"),height=6,width=6)
par(mar=c(0,0,0,0),mfrow=c(3,2),oma=c(2,2,0,0))
setwd(sname)
source('interpolation.R')
setwd(fname)
dev.off()
# ---------------------------------------------------------------------------