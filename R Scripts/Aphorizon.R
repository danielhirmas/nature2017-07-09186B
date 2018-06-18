# ---------------------------------------------------------------------------
# Examine the relationship between climate and ep for A horizons (see
# the selection of A horizons below) that have their mid horizon depth within 
# the upper 25 cm
# ---------------------------------------------------------------------------
setwd(fname)
# Select the Ap horizons
mindepth <- 25 # selects horizons that are less than or equal to this depth 
# (cm)
horizonhold <- c("ABp","Ap","Ap1","AP1","Ap1 / Ap2","Ap1/Ap2","Ap12",
                 "Ap2","Ap3","AP3","Ap4","Apc","Apg","Apt")
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
noAp <- length(hhold)
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
dAp <- d # saves the data frame for this selection
# ---------------------------------------------------------------------------