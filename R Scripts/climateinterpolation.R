# ---------------------------------------------------------------------------
# Interpolation for precip magnitude, precip timing
# D.R. Hirmas
# 10/8/2017
# ---------------------------------------------------------------------------
gd <- al
names(gd)[which(names(gd)=="lat")] <- "y"
names(gd)[which(names(gd)=="long")] <- "x"
#
# Only select the magnitude
gd <- data.frame(x=gd$x,y=gd$y,mag=gd$magnitude)
#
mag.loess <- loess(mag~x*y,data=gd,degree=2,span=0.25,normalize=F)
mag.mar <- list(x = seq(-125,-61,1), y = seq(24,55.5,0.5))
mag.lo <- predict(mag.loess, expand.grid(mag.mar))
usmap <- data.frame(x=map('usa',plot=F)$x[1:6886],
	y=map('usa',plot=F)$y[1:6886])
mag.map <- maskmatrix(mag.lo,usmap)
mag.map[46,9:17] <- NA
#
# Plot
map('usa')
title(main="Precipitation Magnitude",line=0.5)
hsvn <- 100 # number of divisions for the hsv.ramp linear color ramp
hsv.ramp <- hsv(h=rep(235/360,hsvn),s=seq(0, 1, length.out =hsvn),
                v=rep(1,hsvn))
image(mag.mar$x,mag.mar$y,mag.map,col=hsv.ramp,add=T)
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
# Outer margin labels
mtext(expression(,paste("Longitude (",degree,")",sep="")),
side=1,line=0.5,adj=0.49,cex=0.9,outer=TRUE)
mtext(expression(,paste("Latitude (",degree,")",sep="")), side=2,line=0.5,adj=0.49,cex=0.9,outer=TRUE)
#
# Add legend
brk <- c(min(mag.map,na.rm=T),(min(mag.map,na.rm=T)+
max(mag.map,na.rm=T))/2,max(mag.map,na.rm=T))
image.plot(mag.mar$x,mag.mar$y,mag.map,col=hsv.ramp,add=T,
legend.shrink=0.4,legend.only=T,smallplot=c(0.75,0.77,0.24,0.39),
axis.args=list(at=brk,labels=as.character(round(brk,0))))
# smallplot coordinates are left,right,bottom,top in NDC units (0-1)
#
#
# Only select the timing
gd <- al
names(gd)[which(names(gd)=="lat")] <- "y"
names(gd)[which(names(gd)=="long")] <- "x"
gd <- data.frame(x=gd$x,y=gd$y,tim=gd$timing)
#
tim.loess <- loess(tim~x*y,data=gd,degree=2,span=0.25,normalize=F)
tim.mar <- list(x = seq(-125,-61,1), y = seq(24,55.5,0.5))
tim.lo <- predict(tim.loess, expand.grid(tim.mar))
usmap <- data.frame(x=map('usa',plot=F)$x[1:6886],
	y=map('usa',plot=F)$y[1:6886])
tim.map <- maskmatrix(tim.lo,usmap)
tim.map[46,9:17] <- NA
#
# Plot
map('usa')
title(main="Precipitation Timing",line=0.5)
image(tim.mar$x,tim.mar$y,tim.map,col=hsv.ramp,add=T)
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
# Add legend
brk <- c(min(tim.map,na.rm=T),(min(tim.map,na.rm=T)+
max(tim.map,na.rm=T))/2,max(tim.map,na.rm=T))
image.plot(tim.mar$x,tim.mar$y,tim.map,col=hsv.ramp,add=T,
legend.shrink=0.4,legend.only=T,smallplot=c(0.75,0.77,0.24,0.39),
axis.args=list(at=brk,labels=as.character(round(brk,1))))
# smallplot coordinates are left,right,bottom,top in NDC units (0-1)
# ---------------------------------------------------------------------------






