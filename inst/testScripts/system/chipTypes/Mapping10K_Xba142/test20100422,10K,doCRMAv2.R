##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
subset <- c(2,5,3);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doCRMAv2() on an AffymetrixCelSet object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);

# Run CRMAv2 on a subset of the arrays.  This is possible because
# CRMAv2 is truly a single-array method.
ds1 <- doCRMAv2(csR, arrays=subset, verbose=log);
print(ds1);

# Assert that the correct arrays have been processed
stopifnot(getNames(ds1) == getNames(csR)[subset]);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doCRMAv2() by name
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds2 <- doCRMAv2(dataSet, chipType=chipType, arrays=subset, verbose=log);
print(ds2);
stopifnot(getNames(ds2) == getNames(csR)[subset]);

stopifnot(equals(ds2, ds1));
