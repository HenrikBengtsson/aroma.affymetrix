library("aroma.affymetrix");

log <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSetName <- "HapMap270,6.0,CEU,testSet";
chipType <- "GenomeWideSNP_6,Full";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup of annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CDF
cdf <- AffymetrixCdfFile$byChipType(chipType);

# Assert existence of probe-sequence annotation files
acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSetName, cdf=cdf);
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR);
print(acc);
csC <- process(acc, verbose=log);
print(csC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Base-position normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bpn <- BasePositionNormalization(csC, shift=+300);
print(bpn);

csN <- process(bpn, verbose=log);
print(csN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allele-specific chip effect estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE);
print(plm);

if (length(findUnitsTodo(plm)) > 0) {
   # Fit CN probes quickly (~5-10s/array + some overhead)
  units <- fitCnProbes(plm, verbose=log);
  str(units);
  # int [1:945826] 935590 935591 935592 935593 935594 935595 ...

  # Fit remaining units, i.e. SNPs (~5-10min/array)
  units <- fit(plm, verbose=log);
  str(units);
}

ces <- getChipEffectSet(plm);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces);
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract (thetaA, thetaB) for copy-neutral chromosomes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- getCdf(cesN);
gi <- getGenomeInformation(cdf);
units <- getUnitsOnChromosomes(gi, 1:22);
theta <- extractTheta(cesN, units=units);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Estimate copy numbers on the natural scale
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thetaR <- rowMedians(theta[,1,]+theta[,2,], na.rm=TRUE);
C <- 2*(theta[,1,]+theta[,2,])/thetaR;
CA <- 2*theta[,1,]/thetaR;
CB <- 2*theta[,2,]/thetaR;
freqB <- CB/C;

mu <- colMedians(C, na.rm=TRUE);
muA <- colMedians(CA, na.rm=TRUE);
muB <- colMedians(CB, na.rm=TRUE);

Clim <- c(-0.2,3.2);
xlab <- "Freq B (thetaB/theta)";
Clab <- "Copy number";
CAlab <- "Copy number (Allele A)";
CBlab <- "Copy number (Allele B)";

fig <- 1;
if (!devIsOpen(fig <- fig + 1)) {
  devNew();
  layout(matrix(1:9, nrow=3, byrow=TRUE));
  par(mar=c(4,4,0.5,0.5)+0.1);
  centers <- matrix(c(0,2, 1,1, 2,0), nrow=3, ncol=2, byrow=TRUE);
  for (cc in 1:ncol(C)) {
    xx <- CA[,cc];
    yy <- CB[,cc];
    X <- cbind(xx,yy);
    ok <- (is.finite(X) & -1 < X & X < 30);
    ok <- ok[,1] & ok[,2];
    X <- X[ok,];
    smoothScatter(X, xlim=Clim, ylim=Clim, xlab=CAlab, ylab=CBlab);
    abline(h=0:3, lty=3, lwd=1);
    abline(v=0:3, lty=3, lwd=1);
  
    # Plot centers
    fit <- kmeans(X, centers=centers);
    print(fit$centers);
    points(fit$centers, pch=19, col="red");
  }
  devDone();
}

if (!devIsOpen(fig <- fig + 1)) {
  devNew();
  layout(matrix(1:9, nrow=3, byrow=TRUE));
  par(mar=c(4,4,0.5,0.5)+0.1);
  centers <- matrix(c(0,2, 1/2,2, 1,2), nrow=3, ncol=2, byrow=TRUE);
  for (cc in 1:ncol(C)) {
    xx <- freqB[,cc];
    yy <- C[,cc];
    X <- cbind(xx,yy);
    ok <- (is.finite(X) & -1 < X & X < 30);
    ok <- ok[,1] & ok[,2];
    X <- X[ok,];
    smoothScatter(X, xlim=c(0,1), ylim=Clim, xlab=xlab, ylab=Clab);
    abline(h=0:3, lty=3, lwd=1);
    abline(v=0:2/2, lty=3, lwd=1);
  
    # Plot centers
    fit <- kmeans(X, centers=centers);
    print(fit$centers);
    points(fit$centers, pch=19, col="red");
  }
  devDone();
}
