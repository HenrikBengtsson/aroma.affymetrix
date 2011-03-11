###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the CRLMM
# genotype estimates as estimated by oligo.
# It verifies that they give the same results whether or not one
# is normalizing towards the HapMap reference (as defined by oligo).
#
# Author: Henrik Bengtsson
# Created: 2008-12-04
# Last modified: 2009-01-10
#
# Data set:
#  rawData/
#   HapMap270,500K,CEU,testSet/
#     Mapping250K_Nsp/
#       NA06985,Hind,B5,3005533.CEL
#       NA06991,Hind,B6,3005533.CEL
#       NA06993,Hind,B4,4000092.CEL
#       NA06994,Hind,A7,3005533.CEL
#       NA07000,Hind,A8,3005533.CEL
#       NA07019,Hind,A12,4000092.CEL
###########################################################################
library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap270,500K,CEU,testSet";
chipType <- "Mapping250K_Nsp";

cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRLMM according to aroma.affymetrix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- justSNPRMA(csR, normalizeToHapmap=TRUE, returnESet=FALSE, verbose=log);
print(ces);

recalibrate <- TRUE;
crlmm <- CrlmmModel(ces, tags="*,oligo", recalibrate=recalibrate);
print(crlmm);

units <- fit(crlmm, ram="oligo", verbose=log);
str(units);

callSet <- getCallSet(crlmm);
print(callSet);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRLMM according to oligo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("oligoData", getFullName(csR), 
                               getChipType(csR, fullname=FALSE));
path <- Arguments$getWritablePathname(path);
if (!isDirectory(path)) {
  oligo:::justCRLMMv2(getPathnames(csR), tmpdir=path, recalibrate=recalibrate, balance=1.5, verbose=TRUE);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare genotype calls
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
units <- indexOf(cdf, pattern="^SNP");
unitNames <- getUnitNames(cdf, units=units);
units <- units[order(unitNames)];
calls <- extractGenotypes(callSet, units=units, encoding="oligo");
dimnames(calls) <- NULL;

calls0 <- readSummaries("calls", path);
dimnames(calls0) <- NULL;

count <- 0;
for (cc in 1:ncol(calls)) {
  idxs <- whichVector(calls[,cc] != calls0[,cc]);
  count <- count + length(idxs);
  cat(sprintf("%s: ", getNames(callSet)[cc]));
  if (length(idxs) > 0) {
    map <- c("AA", "AB", "BB");
    cat(paste(map[calls[idxs,cc]], map[calls0[idxs,cc]], sep="!="), sep=", ");
  }
  cat("\n");
}
cat(sprintf("Averages number of discrepances per array: %.1f\n", count/ncol(calls)));
errorRate <- count/length(calls);
cat(sprintf("Concordance rate: %.5f%%\n", 100*(1-errorRate)));

stopifnot(errorRate < 1e-4);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare confidence scores
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
confSet <- getConfidenceScoreSet(crlmm);
conf <- extractMatrix(confSet, units=units);
dimnames(conf) <- NULL;

conf0 <- readSummaries("conf", path);
dimnames(conf) <- NULL;

delta <- conf - conf0;
avgDelta <- mean(abs(delta), na.rm=TRUE);
stopifnot(avgDelta < 1e-4);

subplots(ncol(conf));
par(mar=c(3,2,1,1)+0.1);
lim <- c(0,1);
for (cc in 1:ncol(conf)) {
  plot(NA, xlim=lim, ylim=lim);
  abline(a=0, b=1, col="#999999");
  points(conf[,cc], conf0[,cc], pch=".", cex=3);
  rho <- cor(conf[,cc], conf0[,cc]);
  stext(side=3, pos=0, line=-1, sprintf("rho=%.4f", rho));
  stext(side=3, pos=0, getNames(confSet)[cc]);
  cat(sprintf("Array #%d: Correlation: %.4f\n", cc, rho));
}
