# ---------------------------------------------------------------------------
# Interpolation resepx data from R. Kerry
# D.R. Hirmas
# 3/22/2018
# ---------------------------------------------------------------------------
# Save the lowest point in Florida for removal below (i.e., FLremove==T
# is set above)
FLmin <- min(d$lat[which(d$state_admindiv_code=="FL")])

gd <- d
names(gd)[which(names(gd)=="lat")] <- "y"
names(gd)[which(names(gd)=="long")] <- "x"
#
# Select the effective porosity that was calculated using the spatial 
# model
gd <- data.frame(x=gd$x,y=gd$y,resep=gd$RESEPx)
#
resep.loess <- loess(resep~x*y,data=gd,degree=2,span=0.25,normalize=F)
resep.mar <- list(x = seq(-125,-61,1), y = seq(24,55.5,0.5))
resep.lo <- predict(resep.loess, expand.grid(resep.mar))
usmap <- data.frame(x=map('usa',plot=F)$x[1:6886],
	y=map('usa',plot=F)$y[1:6886])
resep.map <- maskmatrix(resep.lo,usmap,FL=FLremove,FLmin=FLmin)
resep.map[46,9:17] <- NA
legresep.map <- adjleg(resep.map,minvalue=minvalueresep,
maxvalue=maxvalueresep)
#
# Plot
map('usa')
title(main="Residual Effective Porosity",line=0.5)
#
# These lines are for an anomaly map centered at zero
hsvn <- 50 # number of divisions for the hsv.ramp linear color ramp
hsv.ramp <- hsv(h=rep(235/360,hsvn),s=seq(0, 1, length.out =hsvn),
                v=rep(1,hsvn)) # blue
hsv.ramp <- rev(c(rev(hsv.ramp),hsv(h=rep(360/360,hsvn),s=seq(0, 100/100, 
                                                          length.out =hsvn),
                v=rep(95/100,hsvn)))) # red

image(resep.mar$x,resep.mar$y,legresep.map,col=hsv.ramp,add=T)
map("state", interior = FALSE, add = TRUE)
map("state", boundary = FALSE, col="black", add = TRUE)
axis(1,at=c(-120,-110,-100,-90,-80,-70))
axis(2,at=c(25,30,35,40,45))
box()
#
# Add points
point <- F # set to false if the addition of points is not desired
if (point==T) {
lhold <- d$long+d$lat
select <- NULL
for (i in 1:length(unique(lhold))) {
	select <- c(select,which(lhold==unique(lhold)[i])[1])
}
points(d$long[select],d$lat[select],pch=1,cex=0.2,col="black")
}
#
# Hide the extra points for the scale
rect(-121,29,-117,31,col="white",lty=0)
#
# Outer margin labels
mtext(expression(,paste("Longitude (",degree,")",sep="")),
side=1,line=0.5,adj=0.49,cex=0.9,outer=TRUE)
mtext(expression(,paste("Latitude (",degree,")",sep="")), side=2,line=0.5,adj=0.49,cex=0.9,outer=TRUE)
#
# Add legend
brk <- c(min(legresep.map,na.rm=T),#(min(legresep.map,na.rm=T)+
#max(legresep.map,na.rm=T))/2,
0,max(legresep.map,na.rm=T))
image.plot(resep.mar$x,resep.mar$y,legresep.map,col=hsv.ramp,add=T,
legend.shrink=0.4,legend.only=T,smallplot=c(0.75,0.77,0.22,0.37),
axis.args=list(at=brk,labels=as.character(round(brk,2))))
# smallplot coordinates are left,right,bottom,top in NDC units (0-1)
# ---------------------------------------------------------------------------






