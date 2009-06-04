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
# Last modified: 2008-12-04
#
# Data set:
#  rawData/
#   HapMap270,100K,CEU,testSet/
#     Mapping50K_Hind240/
#       NA06985,Hind,B5,3005533.CEL
#       NA06991,Hind,B6,3005533.CEL
#       NA06993,Hind,B4,4000092.CEL
#       NA06994,Hind,A7,3005533.CEL
#       NA07000,Hind,A8,3005533.CEL
#       NA07019,Hind,A12,4000092.CEL
###########################################################################

library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

chipType <- "Mapping50K_Hind240";

# Assert that oligo and the correct Platform Design package is installed
library("oligo");
pdPkgName <- cleanPlatformName(chipType);
library(pdPkgName, character.only=TRUE);


normalizeToHapmap <- TRUE;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA according to aroma.affymetrix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName("HapMap270,100K,CEU,testSet", cdf=cdf);
print(csR);

if (normalizeToHapmap) {
  # Load target from PD package
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

qn <- QuantileNormalization(csR, targetDistribution=target, 
                            typesToUpdate="pm", tags=c("*", refTag));
print(qn);
csN <- process(qn, verbose=log);
print(csN);

plm <- RmaSnpPlm(csN, mergeStrands=FALSE, flavor="oligo");
print(plm);
fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);

units <- indexOf(cdf, pattern="^SNP");
unitNames <- getUnitNames(cdf, units=units);
o <- order(unitNames);
units <- units[o];
unitNames <- unitNames[o];
theta <- extractTheta(ces, units=units);
theta <- log2(theta);
dimnames(theta)[[1]] <- unitNames;
rm(o, unitNames);

# Make sure pairs are order as (sense, antisense)
dirs <- getGroupDirections(cdf, units=units);
dirs <- matrix(unlist(dirs, use.names=FALSE), nrow=4);
idxs <- which(dirs[1,] == 2);
theta[idxs,,] <- theta[idxs,c(3,4,1,2),];
rm(units,dirs,idxs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA according to oligo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
eSet <- justSNPRMA(getPathnames(csR), normalizeToHapmap=normalizeToHapmap, verbose=TRUE);
print(eSet);



# Extract theta array
naValue <- as.double(NA);
theta0 <- array(naValue, dim=c(nrow(eSet), 4, ncol(eSet)));
theta0[,1,] <- senseThetaA(eSet);
theta0[,2,] <- senseThetaB(eSet);
theta0[,3,] <- antisenseThetaA(eSet);
theta0[,4,] <- antisenseThetaB(eSet);
dimnames(theta0) <- list(NULL, NULL, NULL);
dimnames(theta0)[[1]] <- featureNames(eSet);
dimnames(theta0)[[3]] <- getNames(csR);


# Assert that the dimensions are the same
stopifnot(identical(dim(theta), dim(theta0)));

# Assert that the ordering of units and arrays are the same
stopifnot(identical(dimnames(theta), dimnames(theta0)));

# Assert that the estimates are very similar
tol <- 1e-4;
res <- all.equal(theta, theta0, tolerance=tol);
print(res);
stopifnot(res);
