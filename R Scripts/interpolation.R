# ---------------------------------------------------------------------------
# Interpolation for EP
# D.R. Hirmas
# 10/8/2017
# ---------------------------------------------------------------------------
gd <- d
names(gd)[which(names(gd)=="lat")] <- "y"
names(gd)[which(names(gd)=="long")] <- "x"
#
# Only select the effective porosity
gd <- data.frame(x=gd$x,y=gd$y,ep=gd$ep)
#
ep.loess <- loess(ep~x*y,data=gd,degree=2,span=0.25,normalize=F)
ep.mar <- list(x = seq(-125,-61,1), y = seq(24,55.5,0.5))
ep.lo <- predict(ep.loess, expand.grid(ep.mar))
usmap <- data.frame(x=map('usa',plot=F)$x[1:6886],
	y=map('usa',plot=F)$y[1:6886])
ep.map <- maskmatrix(ep.lo,usmap,FL=FLremove,FLmin=FLmin)
ep.map[46,9:17] <- NA
legep.map <- adjleg(ep.map,minvalue=minvalueep,maxvalue=maxvalueep)
#
# Plot
map('usa')
title(main="Effective Porosity",line=0.5)
hsvn <- 100 # number of divisions for the hsv.ramp linear color ramp
# Use this line of code for a red linear color ramp:
hsv.ramp <- hsv(h=rep(130/360,hsvn),s=seq(0, 90/100, length.out =hsvn),
                v=rep(45/100,hsvn),alpha=seq(0,1,length.out=hsvn)) # green
image(ep.mar$x,ep.mar$y,legep.map,col=hsv.ramp,add=T)
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
brk <- c(min(legep.map,na.rm=T),(min(legep.map,na.rm=T)+
max(legep.map,na.rm=T))/2,max(legep.map,na.rm=T))
image.plot(ep.mar$x,ep.mar$y,legep.map,col=hsv.ramp,add=T,
legend.shrink=0.4,legend.only=T,smallplot=c(0.75,0.77,0.22,0.37),
axis.args=list(at=brk,labels=as.character(round(brk,2))))
# smallplot coordinates are left,right,bottom,top in NDC units (0-1)
# ---------------------------------------------------------------------------






