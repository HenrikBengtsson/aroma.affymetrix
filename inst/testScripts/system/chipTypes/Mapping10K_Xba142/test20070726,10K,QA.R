library("aroma.affymetrix")
log <- Verbose(threshold=-4, timestamp=TRUE);


dataSetName <- "Jeremy_2007-10k";
chipType <- "Mapping10K_Xba142";

# Expected sample names
sampleNames <- c("0001-7", "0002-10", "0004-13", "0005-14", "0007-18", 
                      "0008-19", "0010-22", "2-DPrrr", "MH12", "MH18");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
keep <- 1:6;
csR <- extract(csR, keep);
sampleNames <- sampleNames[keep];
print(csR);
stopifnot(identical(getNames(csR), sampleNames));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fitting log-additive model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaPlm(csR);
print(plm);
fit(plm, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Some basic quality scores
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- getChipEffectSet(plm);

# Boxplots of log2(theta), RLE, and NUSE
layout(matrix(1:4, ncol=2, byrow=TRUE));
plotBoxplot(ces, type="theta", transform=log2);
plotBoxplot(ces, type="RLE", arrays=c(2,4:6));
plotBoxplot(ces, type="NUSE");
devDone();


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Advanced usage
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculating statistics without plotting
theta <- boxplotStats(ces, type="theta", transform=log2);
nuse <- boxplotStats(ces, type="NUSE");
rle <- boxplotStats(ces, type="RLE");

# Subset of arrays (avoids calculating stats for all arrays)
rleB <- boxplotStats(ces, type="RLE", arrays=c(2,4:6));
for (name in names(rleB)) {
  stopifnot(identical(rleB[name], rle[name]));
}

# Plotting the above statistics
layout(matrix(1:4, ncol=2, byrow=TRUE));
plotBoxplotStats(theta, main="theta");
plotBoxplotStats(rle[c(2,4:6)], main="RLE");
plotBoxplotStats(nuse, main="NUSE");
devDone();


# Calculates unit-specific RLE and NUSE scores
units <- 1000+1:5000;
theta <- extractMatrix(ces, units=units);
rle <- extractMatrix(ces, units=units, field="RLE");
nuse <- extractMatrix(ces, units=units, field="NUSE");

# Plotting the above statistics
layout(matrix(1:4, ncol=2, byrow=TRUE));
plotDensity(log2(theta), main="theta");
plotDensity(rle[,c(2,4:6)], main="RLE");
plotDensity(nuse, main="NUSE");
devDone();

# ...same, but basic unit annotation data added
units <- 1000+1:500;
theta <- extractDataFrame(ces, units=units, addNames=TRUE);
rle <- extractDataFrame(ces, units=units, field="RLE", addNames=TRUE);
nuse <- extractDataFrame(ces, units=units, field="NUSE", addNames=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Robust estimate of standard deviation of raw CNs
# (The default is to use a first-order difference variance estimator)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (FALSE) {
  # Estimate it for some autosomal chromosomes and ChrX
  # From aroma.core v1.3.3 it is no longer possible to define
  # a segmentation model if the data set is not a copy number
  # data set. /HB 2009-11-19
  cns <- CbsModel(ces);
  res <- estimateSds(cns, chromosomes=c(1:6, 23), verbose=log);
  chrX <- which(rownames(res) == "23");
  
  layout(matrix(1:2, nrow=2));
  par(mar=c(4,4,0.5,2)+0.1);
  
  xlim <- c(1,(ncol(res)+2));
  ylim <- c(0, 1.05*max(res, na.rm=TRUE));
  ylab <- expression(hat(sigma)==s[Delta](log2(theta/theta[R])));
  cols <- seq(length=nrow(res));
  ltys <- rep(1, times=nrow(res));
  cols[chrX] <- "blue";
  ltys[chrX] <- 4;
  
  plot(NA, xlim=xlim, ylim=ylim, xlab="", ylab=ylab);
  for (kk in seq(length=nrow(res))) {
    chromosome <- rownames(res)[kk];
    points(res[kk,], col=cols[kk], pch=19);
    lines(res[kk,], col=cols[kk], lty=ltys[kk], lwd=2);
  }
  legend("topright", col=cols, lty=ltys, pch=19, lwd=2, 
                     legend=sprintf("Chr%s", rownames(res)), bty="n");
  
  
  df <- as.data.frame(res[-chrX,,drop=FALSE]);
  colnames(df) <- seq(length=nrow(df));
  boxplot(df, ylim=ylim, ylab=ylab, xlab="Array");
  points(res[chrX,], col="blue", pch=19, cex=1.5);
  legend("bottomright", col=c("black", cols[chrX]), pch=19, lwd=2, 
                     legend=c("Autosomal", "ChrX"), horiz=TRUE, bty="n");
  devDone();
} # if (FALSE)
