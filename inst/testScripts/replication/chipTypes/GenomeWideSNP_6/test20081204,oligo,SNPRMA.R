###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the SNPRMA 
# chip-effect estimates as estimated by oligo.
# It verifies that they give the same results whether or not one
# is normalizing towards the HapMap reference (as defined by oligo).
#
# Author: Henrik Bengtsson
# Created: 2008-12-04
# Last modified: 2010-05-09
#
# Data set:
#  rawData/
#   HapMap270,6.0,CEU,testSet/
#     GenomeWideSNP_6/
#       NA06985.CEL
#       NA06991.CEL
#       NA06993.CEL
#       NA06994.CEL
#       NA07000.CEL
#       NA07019.CEL
###########################################################################
library("aroma.affymetrix");
library("oligo");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

chipType <- "GenomeWideSNP_6";

normalizeToHapmap <- TRUE;
normalizeToHapmap <- FALSE;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pkg <- Package("oligo");
if (isOlderThan(pkg, "1.8.0")) {
  # oligo v1.7.x and older
  thetaA <- oligo:::thetaA;
  thetaB <- oligo:::thetaB;
} else if (isOlderThan(pkg, "1.12.0")) {
  # oligo v1.11.x and older
  # Nothing to do
} else {
  # oligo v1.12.0 and newer
  thetaA <- function(eSet) {
    ad <- assayData(eSet);
    ad$alleleA;
  }
  thetaB <- function(eSet) {
    ad <- assayData(eSet);
    ad$alleleB;
  }
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA according to aroma.affymetrix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName("HapMap270,6.0,CEU,testSet", cdf=cdf);
print(csR);

if (normalizeToHapmap) {
  # Load target from PD package
  pdPkgName <- oligo::cleanPlatformName(chipType);
  stopifnot(isPackageInstalled(pdPkgName));
  path <- system.file("extdata", package=pdPkgName);
  filename <- sprintf("%sRef.rda", pdPkgName);
  pathname <- Arguments$getReadablePathname(filename, path=path, 
                                                    mustExits=TRUE);
  target <- loadToEnv(pathname)$reference;
  refTag <- "HapMapRef";
} else {
  target <- NULL;
  refTag <- NULL;
}

# justSNPRMA() operates only on *SNP* units (CN units ignored).
# For this reason we here *estimate* the normalization function based
# on these units only, but for convenience we will apply it to all
# units including CN units.
units <- indexOf(cdf, pattern="^SNP");
cells <- getCellIndices(cdf, units=units, unlist=TRUE, useNames=FALSE);
qn <- QuantileNormalization(csR, targetDistribution=target, 
                                subsetToAvg=cells, typesToUpdate="pm", 
                                         tags=c("*", "SNPs", refTag));
print(qn);
rm(cells);
csN <- process(qn, verbose=log);
print(csN);

plm <- RmaSnpPlm(csN, mergeStrands=FALSE, flavor="oligo");
print(plm);
fitCnProbes(plm, verbose=log);
fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);

unitNames <- getUnitNames(cdf, units=units);
o <- order(unitNames);
unitNames <- unitNames[o];
units <- units[o];
theta <- extractTheta(ces, groups=1:2, units=units, verbose=log);
theta <- log2(theta);
dimnames(theta)[[1]] <- unitNames;
rm(o, units, unitNames);


# For debugging purposes
print(sessionInfo());

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA according to oligo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
eSet <- justSNPRMA(getPathnames(csR), normalizeToHapmap=normalizeToHapmap, verbose=TRUE);
print(eSet);

# Extract theta array
naValue <- as.double(NA);
theta0 <- array(naValue, dim=c(nrow(eSet), 2, ncol(eSet)));
theta0[,1,] <- thetaA(eSet);
theta0[,2,] <- thetaB(eSet);
dimnames(theta0) <- list(NULL, NULL, NULL);
dimnames(theta0)[[1]] <- featureNames(eSet);
dimnames(theta0)[[3]] <- getNames(csR);

# Assert that the dimensions are the same
stopifnot(identical(dim(theta), dim(theta0)));

# Assert that the ordering of units and arrays are the same
stopifnot(identical(dimnames(theta), dimnames(theta0)));

# Assert that the estimates are very similar
tol <- 0.012;
stopifnot(all.equal(theta, theta0, tolerance=tol));
