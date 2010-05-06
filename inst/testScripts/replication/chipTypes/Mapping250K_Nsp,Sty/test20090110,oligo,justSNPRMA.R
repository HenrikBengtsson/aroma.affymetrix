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
# Last modified: 2010-05-06
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

# Local functions
compareESets <- function(eSet1, eSet2, FUN=NULL, tolerance=1e-4) {
  pkg <- Package("oligo");
  if (!isOlderThan(pkg, "1.12.0")) {
    # For oligo v1.12.0 or newer
    stopifnot(all.equal(featureNames(eSet1), featureNames(eSet2)));

    ad1 <- assayData(eSet1);
    l1 <- as.list(ad1);
    ad2 <- assayData(eSet2);
    l2 <- as.list(ad2);
    stopifnot(all.equal(l1, l2, tolerance=tolerance));
  } else {
    # For oligo v1.11.x or older
    if (is.null(FUN)) {
      fcns <- list(senseThetaA, senseThetaB, 
                   antisenseThetaA, antisenseThetaB, 
                   featureNames);
      for (fcn in fcns) {
        compareESets(eSet1, eSet2, FUN=fcn, tolerance=tolerance);
      }
    } else {
      stopifnot(all.equal(FUN(eSet1), FUN(eSet2), tolerance=tolerance));
    }
  }
  TRUE;
} # compareESets()

dataSet <- "HapMap270,500K,CEU,testSet";
chipType <- "Mapping250K_Nsp";

# Assert that oligo and the correct Platform Design package is installed
library("oligo");
pdPkgName <- cleanPlatformName(chipType);
library(pdPkgName, character.only=TRUE);


normalizeToHapmap <- TRUE;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA according to aroma.affymetrix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);

eSet <- justSNPRMA(csR, normalizeToHapmap=normalizeToHapmap, verbose=log);
print(eSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA according to oligo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
eSet0 <- justSNPRMA(getPathnames(csR), normalizeToHapmap=normalizeToHapmap, verbose=TRUE);
print(eSet0);

stopifnot(compareESets(eSet, eSet0));
