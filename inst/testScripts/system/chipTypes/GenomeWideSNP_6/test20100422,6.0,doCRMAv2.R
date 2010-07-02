##########################################################################
# Data set:
# ...
##########################################################################
library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "HapMap270,6.0,CEU,testSet";
chipType <- "GenomeWideSNP_6,Full";
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
# doCRMAv2() by name and CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds2 <- doCRMAv2(dataSet, cdf=cdf, arrays=subset, verbose=log);
print(ds2);
stopifnot(getNames(ds2) == getNames(csR)[subset]);
stopifnot(equals(ds2, ds1));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doCRMAv2() by name and chip type (with chiptype tags)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds3 <- doCRMAv2(dataSet, chipType=chipType, arrays=subset, verbose=log);
print(ds3);
stopifnot(getNames(ds3) == getNames(csR)[subset]);
stopifnot(equals(ds3, ds1));
