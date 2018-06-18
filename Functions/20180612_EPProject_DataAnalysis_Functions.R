# ---------------------------------------------------------------------------
# Functions for effective porosity manuscript
# D.R. Hirmas
# ---------------------------------------------------------------------------
#
# ---------------------------------------------------------------------------
# Function to mask a matrix outside the borders of a polygon
# ---------------------------------------------------------------------------
# Change the working directory
maskmatrix <- function(mat,gdf,FL=F,FLmin=25) {
	# mat is the matrix created from loess.predict fit component to be masked;
	# gdf is a data.frame with x and y values representing the vertices of 
	# a polygon
	# FL is logical to indicate whether to remove (T) the interpolation
	# below the lowest recorded data point whose latitude is given by FLmin
	dn <- dimnames(mat)
	xcoord <- dn$x
	x <- NULL
	for (i in 1:length(xcoord)) {
		x <- c(x,as.numeric(strsplit(xcoord,"x=")[[i]][2]))
	}
	ycoord <- dn$y
	y <- NULL
	for (i in 1:length(ycoord)) {
		y <- c(y,as.numeric(strsplit(ycoord,"y=")[[i]][2]))
	}
	nm <- matrix(0,nrow=length(x),ncol=length(y))
	for (i in 1:length(gdf$x)) {
		xsav <- round(gdf$x[i],0)
		ysav <- round(gdf$y[i],0)
		xc <- which(x == xsav)
		yc <- which(y == ysav)
		nm[xc,yc] <- 1
	}
	hold <- 0
	for (i in 1:nrow(nm)) {
		for (j in 1:ncol(nm)) {
			if (nm[i,j] == 1) {
				if (hold == 0) {
					hold <- j
				} else {
					nm[i,hold:j] <- 1
				}
			}
		}
		hold <- 0
	}
	for (i in 1:nrow(nm)) {
		for (j in 1:ncol(nm)) {
			if (nm[i,j] == 0) mat[i,j] <- NA
		}
	}
	# Remove the Florida extrapolation for the non existing points
	if (FL==T) {
		xhold <- which((x<=-78) & (x>=-85))
		yhold <- which(y<FLmin)
		mat[xhold,yhold] <- NA
	}
	mat
}
# ---------------------------------------------------------------------------
#
#
#
# ---------------------------------------------------------------------------
# Function to add help adjust the scale of the map legends
# ---------------------------------------------------------------------------
# Change the working directory
adjleg <- function(mat,minvalue=NA,maxvalue=NA) {
	# mat is the matrix created from loess.predict fit component and
	# modified by the maskmatrix command
	# minvalue and maxvalue are the values to be assigned to an empty spot 
	# on the map to change the scale
	dn <- dimnames(mat)
	xcoord <- dn$x
	x <- NULL
	for (i in 1:length(xcoord)) {
		x <- c(x,as.numeric(strsplit(xcoord,"x=")[[i]][2]))
	}
	ycoord <- dn$y
	y <- NULL
	for (i in 1:length(ycoord)) {
		y <- c(y,as.numeric(strsplit(ycoord,"y=")[[i]][2]))
	}
	if (is.na(minvalue)==F) { # Note this has to match exactly
		# need to rewrite to allow for the case when the interpolation
		# coordinates are not integers
		xhold <- which(x==-120) 
		yhold <- which(y==30)
		mat[xhold,yhold] <- minvalue
	}
	if (is.na(maxvalue)==F) {
		xhold <- which(x==-118)
		yhold <- which(y==30)
		mat[xhold,yhold] <- maxvalue
	}
	mat
}
# ---------------------------------------------------------------------------
#
#
#
# ---------------------------------------------------------------------------
# Function to predict EEMT from MAT and MAP (Rasmussen and Tabor, 2007)
# ---------------------------------------------------------------------------
EEMT <- function(mat,map) {
  # mat = mean annual temperature (C)
  # map = mean annual precipitation (mm)
  347134*exp(-0.5*(((mat-21.5)/-10.1)^2+((map-4412)/1704)^2)) # EEMT
}
# ---------------------------------------------------------------------------
#
#
#
# ---------------------------------------------------------------------------
# Function to predict Ks from EP, FC, and WP following (Rawls et al., 1998)
# or a similar approach suggested by Han et al. (2008).
# ---------------------------------------------------------------------------
Ksatpred <- function(ep,fc,wp,tp=NA,unit="mmh",method="Rawls") {
  # Uses either the Rawls et al. (1998) or the Han et al. (2008) method
  # For Han's need to specify the total porosity (tp)
  # ep = effective porosity
  # fc = field capacity
  # wp = wilting point
  # Possible values for unit are "mmh" = mm/h, "ums" = um/s,
  # "cmd" = cm/d, "ms" = m/s
  # All porosity and capacity values should be given on a volumetric basis
  if (method=="Rawls") {
    lambda <- abs(log10(wp/fc)/log10(-1500/-33)) # lambda calculated as
    # the slope of the wrc between wp and fc
    D <- 3-lambda
    Ksat <- 1930*ep^D # mm/h
  }
  if (method=="Han") {
    lambda <- abs(log10(tp/fc)/log10(-0.1/-33)) # lambda calculated as
    # the slope of the wrc between tp and fc as suggested by Han et al (2008)
    # Note: this does not appear to be Han's model exactly, but it does use
    # lambda calculated from the same domain as effective porosity which
    # is a similar approach to Han. Also note that the water content at
    # corresponding to the tp point is assigned 0.1 kPa since the log
    # of 0 is undefined.
    D <- 3-lambda
    Ksat <- 1930*ep^D # mm/h
  }
  if (unit=="ums") Ksat <- Ksat*1000/3600 # um/s
  if (unit=="cmd") Ksat <- Ksat/10*24 # cm/d
  if (unit=="ms") Ksat <- Ksat/1000/3600 # m/s
  Ksat
}
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Functions to rotate a matrix
# ---------------------------------------------------------------------------
rotate <- function(x) t(apply(x, 2, rev))
counterclockwise <- function(x) rotate(rotate(rotate(x)))
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Function to plot a filled contour surface on a ternary plot with 
# sand, silt, clay, and z (any fourth variable) data
# ---------------------------------------------------------------------------
txttrifilled <- function(sand,clay,z) {
  # textural data should be provided in percent
  nclay <- clay/100*sqrt(3)/2
  nsand <- 1-nclay/tan(pi/3)-sand/100
  topo <- data.frame(nsand,nclay,z)
# 
  require(matlab)
  names(topo) <- c("x","y","z")
#
  ml <- loess(z ~ x*y, topo)
#
  trian <- expand.grid(base=seq(0,1,l=100*2), high=seq(0,sin(pi/3),l=87*2))
  trian <- subset(trian, (base*sin(pi/3)*2)>high)
  trian <- subset(trian, ((1-base)*sin(pi/3)*2)>high)
#
  new2 <- data.frame(x=trian$base,y=trian$high)
#
  trian$yhat <- predict(ml, newdata=new2)
#
  grade.trellis <- function(from=0.2, to=0.8, step=0.2, col=1, lty=2, 
                            lwd=0.5){
    x1 <- seq(from, to, step)
    x2 <- x1/2
    y2 <- x1*sqrt(3)/2
    x3 <- (1-x1)*0.5+x1
    y3 <- sqrt(3)/2-x1*sqrt(3)/2
    panel.segments(x1, 0, x2, y2, col=col, lty=lty, lwd=lwd)
    panel.text(x1, 0, label=rev(x1)*100, pos=1)
    panel.segments(x1, 0, x3, y3, col=col, lty=lty, lwd=lwd)
    panel.text(x2, y2, label=x1*100, pos=2)
    panel.segments(x2, y2, 1-x2, y2, col=col, lty=lty, lwd=lwd)
    panel.text(x3, y3, label=x1*100, pos=4)
  }
#
  levelplot(yhat~base*high, trian, aspect="iso", cuts=30,
            xlim=c(-0.1,1.1), ylim=c(-0.1,0.96),
            xlab=NULL, ylab=NULL, contour=TRUE,
            par.settings=list(axis.line=list(col=NA),
                              axis.text=list(col=NA),
                              regions=list(col=jet.colors(100))),
            panel=function(..., at, contour=TRUE, labels=NULL){
            panel.levelplot(..., at=at, contour=contour, #labels=labels,
                            lty=3, lwd=1, col="black")
          })
  trellis.focus("panel", 1, 1, highlight=FALSE)
  lpoints(topo$x,topo$y,col="black")
  lpolygon(c(0,0.5,1,0),c(0,sqrt(3)/2,0,0))
  panel.segments(c(0,0,0.5), c(0,0,sqrt(3)/2), c(1,1/2,1), 
                 c(0,sqrt(3)/2,0),lwd=1.5)
  grade.trellis()
  panel.text(1/4, sqrt(3)/4+0.075, label="Clay (%)",pos=2,srt=60,
             offset=1.5,cex=1.2)
  panel.text(3/4, sqrt(3)/4+0.075, label="Silt (%)",pos=4,srt=-60,
             offset=1.5,cex=1.2)
  panel.text(0.5, 0, label="Sand (%)", pos=1,offset=1.7,cex=1.2)
  trellis.unfocus()
}
# ---------------------------------------------------------------------------
#
#
#
#
# ---------------------------------------------------------------------------
# Functions to identify multivariate outliers
# ---------------------------------------------------------------------------
mvOutlier <- function (data, qqplot = TRUE, 
                       method = c("quan", "adj.quan")) { 
  # From the MVN package
  require(robustbase)
  dataframe = as.data.frame(data)
  dname <- deparse(substitute(data))
  method <- match.arg(method)
  n <- dim(data)[1]
  p <- dim(data)[2]
  covr <- covMcd(data, alpha = 0.5)
  mah <- mahalanobis(data, center = covr$center, cov = covr$cov)
  d <- mah
  sortMah <- data.frame(sort(mah, decreasing = TRUE))
  out <- cbind(round(sortMah, 3), NA)
  colnames(out) <- c("MD", "Outlier")
  if (method == "adj.quan") {
    crt <- arw(x = data, m0 = covr$center, c0 = covr$cov, 
               alpha = 0.025)$cn
    for (i in 1:n) {
      {
        if (sortMah[i, ] > crt) {
          out[i, 2] <- "TRUE"
        }
        else {
          out[i, 2] <- "FALSE"
        }
      }
    }
    if (qqplot) {
      d <- mah
      r <- rank(d)
      chi2q <- qchisq((r - 0.5)/n, p)
      colors = NULL
      for (i in 1:n) {
        if (d[i] > crt) 
          colors[i] = "red"
        else colors[i] = "black"
      }
      plot(d, chi2q, pch = 16, main = "Adjusted Chi-Square Q-Q Plot", 
           xlab = "Robust Squared Mahalanobis Distance", 
           ylab = "Chi-Square Quantile", col = colors)
      abline(v = crt, lwd = 2, col = "blue")
      tbl = table(out[, 2])
      legend("topleft", legend = 
               c(paste("Outliers (n=", if (is.na(tbl[2])) 0 else tbl[2], 
                       ")", sep = ""),
                 paste("Non-outliers (n=", if (is.na(tbl[1])) 0 else tbl[1],
                       ")", sep = "")), col = c("red", "black"), 
             pch = 16, bty = "n", )
      if (max(d) >= crt) {
        text(crt - 0.2, 2, paste("Quantile: ", round(crt, 3)), 
             srt = 90, 
             pos = 3, col = "blue")
      }
    }
    newData <- out[out$Outlier %in% "FALSE", ]
    ind <- sort(row.names(newData))
    newData <- data[ind, ]
    result <- list(out, newData)
    names(result) <- c("outlier", "newData")
  }
  if (method == "quan") {
    chiSq <- qchisq(0.975, p)
    for (i in 1:n) {
      {
        if (sortMah[i, ] > chiSq) {
          out[i, 2] <- "TRUE"
        }
        else {
          out[i, 2] <- "FALSE"
        }
      }
    }
    if (qqplot) {
      d <- mah
      r <- rank(d)
      chi2q <- qchisq((r - 0.5)/n, p)
      colors = NULL
      for (i in 1:n) {
        if (d[i] > chiSq) 
          colors[i] = "red"
        else colors[i] = "black"
      }
      plot(d, chi2q, pch = 16, col = colors, main = "Chi-Square Q-Q Plot", 
           xlab = "Robust Squared Mahalanobis Distance", 
           ylab = "Chi-Square Quantile")
      abline(v = chiSq, lwd = 2, col = "red")
      tbl = table(out[, 2])
      legend("topleft", 
             legend = c(paste("Outliers (n=", 
                              if (is.na(tbl[2])) 0 else tbl[2], ")", 
                              sep = ""),
                        paste("Non-outliers (n=", 
                              if (is.na(tbl[1])) 0 else tbl[1],
                              ")", sep = "")), 
             col = c("red", "black"), pch = 16, 
             bty = "n", )
      if (max(d) >= chiSq) {
        text(chiSq - 0.2, 2, paste("Quantile: ", 
                                   round(chiSq,3)), 
             srt = 90, pos = 3, col = "red")
      }
    }
    newData <- out[out$Outlier %in% "FALSE", ]
    ind <- sort(row.names(newData))
    newData <- data[ind, ]
    result <- list(out, newData)
    names(result) <- c("outlier", "newData")
  }
  return(result)
}
#
arw <- function (x, m0, c0, alpha, pcrit) # from the mvoutlier package
{
  n <- nrow(x)
  p <- ncol(x)
  if (missing(pcrit)) {
    if (p <= 10) 
      pcrit <- (0.24 - 0.003 * p)/sqrt(n)
    if (p > 10) 
      pcrit <- (0.252 - 0.0018 * p)/sqrt(n)
  }
  if (missing(alpha)) 
    delta <- qchisq(0.975, p)
  else delta <- qchisq(1 - alpha, p)
  d2 <- mahalanobis(x, m0, c0)
  d2ord <- sort(d2)
  dif <- pchisq(d2ord, p) - (0.5:n)/n
  i <- (d2ord >= delta) & (dif > 0)
  if (sum(i) == 0) 
    alfan <- 0
  else alfan <- max(dif[i])
  if (alfan < pcrit) 
    alfan <- 0
  if (alfan > 0) 
    cn <- max(d2ord[n - ceiling(n * alfan)], delta)
  else cn <- Inf
  w <- d2 < cn
  if (sum(w) == 0) {
    m <- m0
    c <- c0
  }
  else {
    m <- apply(x[w, ], 2, mean)
    c1 <- as.matrix(x - rep(1, n) %*% t(m))
    c <- (t(c1 * w) %*% c1)/sum(w)
  }
  list(m = m, c = c, cn = cn, w = w)
}
#
outlierID <- function(df) {
  # Creates a logical vector to test to see if the point is a bivariate
  # outlier; vector corresponds to the order of the rows in the input df
  # df = data.frame; only the first two columns will be used to detect the
  # outliers with; those first two columns need to be numeric
  #
  # Please note: there is an issue with mvOutlier where if the 
  # row.names of the dataframe are NULL from 1 to nrow(df), the algorithm
  # will designate the number of outliers correctly but will remove the
  # wrong ones in the resulting dfout$newData dataframe.
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
# ---------------------------------------------------------------------------