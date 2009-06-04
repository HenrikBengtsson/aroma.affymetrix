library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

chipType <- "GenomeWideSNP_6,Full";

cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName("HapMap270,6.0,CEU,testSet", cdf=cdf);
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA with no need for PD packages
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
eSet <- justSNPRMA(csR, normalizeToHapmap=FALSE, verbose=log);
print(eSet);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA with need for PD packages
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert that oligo and the correct Platform Design package is installed
require("oligo") || throw("Package not loaded: oligo");
pdPkgName <- cleanPlatformName(getChipType(cdf, fullname=FALSE));
library(pdPkgName, character.only=TRUE);

eSet <- justSNPRMA(csR, normalizeToHapmap=TRUE, verbose=log);
print(eSet);
