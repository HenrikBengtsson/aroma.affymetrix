library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE13372,testset";
chipType <- "GenomeWideSNP_6,Full";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract (CN, freqB)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesN <- res$cesN;
ceR <- getAverage(cesN, verbose=verbose);
print(ceR);

ugp <- getAromaUgpFile(cesN);
print(ugp);

chr <- 2
units <- getUnitsOnChromosome(ugp, chromosome=2, region=c(75,90)*1e6);
pos <- getPositions(ugp, units=units) / 1e6;

thetaR <- extractTotalAndFreqB(ceR, units=units)[,"total"];
ce <- getFile(cesN, 1);
data <- extractTotalAndFreqB(ce, units=units);
data[,"total"] <- 2*data[,"total"] / thetaR;

C <- data[,"total"];
B <- data[,"freqB"];

okC <- whichVector(is.finite(C));
okB <- whichVector(is.finite(B));
fit <- kmeans(B[okB], centers=c(0,1/2,1))

layout(matrix(1:2, ncol=1))
par(mar=c(3,4,2,1)+0.1, pch=".")
xlim <- range(pos, na.rm=TRUE);
plot(NA, xlim=xlim, ylim=c(0,4), xlab="pos", ylab=expression(C))
abline(h=0:4, lty=3, lwd=2, col="#999999")
abline(h=median(C, na.rm=TRUE), lwd=2, col="red")
points(pos,C, cex=3);
lines(smooth.spline(pos,C), lwd=2, col="blue");
abline(v=c(83.1,83.7))
stext(side=3, pos=0, getName(ce))
stext(side=3, pos=1, sprintf("Chr %d", chr))
plot(NA, xlim=xlim, ylim=c(0,1), xlab="pos", ylab=expression(beta))
abline(h=fit$centers, lwd=2, col="black")
points(pos[okB], B[okB], cex=3, col=fit$cluster+1)
abline(v=c(83.1,83.7))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Stratify by SNPs and CN units
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
units <- getUnitsOnChromosome(ugp, chromosome=2);
pos <- getPositions(ugp, units=units) / 1e6;

thetaR <- extractTotalAndFreqB(ceR, units=units)[,"total"];
ce <- getFile(cesN, 1);
data <- extractTotalAndFreqB(ce, units=units);
data[,"total"] <- 2*data[,"total"] / thetaR;

C <- data[,"total"];
B <- data[,"freqB"];

okC <- whichVector(is.finite(C));
okB <- whichVector(is.finite(B));
fit <- kmeans(B[okB], centers=c(0,1/2,1))

isCN <- (getUnitTypes(cdf, units=units) == 5);
layout(matrix(1:3, ncol=1))
par(mar=c(3,4,2,1)+0.1, pch=".")
xlim <- range(pos, na.rm=TRUE);
plot(NA, xlim=xlim, ylim=c(0,4), xlab="pos", ylab=expression(C))
abline(h=0:4, lty=3, lwd=2, col="#999999")
abline(h=median(C, na.rm=TRUE), lwd=2, col="red")
x <- pos[!isCN]; y <- C[!isCN];
points(x,y, cex=3);
lines(smooth.spline(x,y), lwd=2, col="blue");
abline(v=c(83.1,83.7))
stext(side=3, pos=0, getName(ce))
stext(side=3, pos=1, sprintf("Chr %d", chr))
plot(NA, xlim=xlim, ylim=c(0,4), xlab="pos", ylab=expression(C))
abline(h=0:4, lty=3, lwd=2, col="#999999")
abline(h=median(C, na.rm=TRUE), lwd=2, col="red")
x <- pos[isCN]; y <- C[isCN];
points(x,y, cex=3);
lines(smooth.spline(x,y), lwd=2, col="blue");
abline(v=c(83.1,83.7))
stext(side=3, pos=0, getName(ce))
stext(side=3, pos=1, sprintf("Chr %d", chr))
plot(NA, xlim=xlim, ylim=c(0,1), xlab="pos", ylab=expression(beta))
abline(h=fit$centers, lwd=2, col="black")
points(pos[okB], B[okB], cex=3, col=fit$cluster+1)
abline(v=c(83.1,83.7))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot (lambda, C)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
isCN <- (getUnitTypes(cdf, units=units) == 5);
fl <- getFragmentLengths(si, units=units);

xlim <- c(0,2000);
Clim <- c(0,4);
layout(matrix(1:4, ncol=1));
par(mar=c(1,1,1,1)+0.1, pch=".")
for (kk in 1:2) {
  if (kk == 1) {
    idxs <- whichVector(isCN);
  } else {
    idxs <- whichVector(!isCN);
  }
  for (ee in 1:2) {
    x <- fl[idxs,ee];
    y <- C[idxs];
    ok <- whichVector(is.finite(x) & is.finite(y));
    x <- x[ok];
    y <- y[ok];
    plot(x,y, cex=2, xlim=xlim, ylim=Clim);
    lines(smooth.spline(x,y), lwd=2, col="red");
  }
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot density of freqB stratified by CN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
breaks <- quantile(C, probs=seq(0,1,by=0.1), na.rm=TRUE)
cuts <- cut(C, breaks=c(0,1.4,1.68,10));
cuts <- as.integer(cuts);
ucuts <- unique(cuts);
cols <- terrain.colors(length(ucuts));
plot(NA, xlim=c(0,1), ylim=c(0,3), xlab="B", ylab="density");
for (kk in seq(along=ucuts)) {
  idxs <- whichVector(cuts == kk);
  zkk <- na.omit(B[idxs]);
  if (length(zkk) > 1) {
    print(length(zkk));
    d <- density(zkk, adjust=0.6, kernel="biweight");
    lines(d, lwd=2, col=cols[kk]);
  }
}
