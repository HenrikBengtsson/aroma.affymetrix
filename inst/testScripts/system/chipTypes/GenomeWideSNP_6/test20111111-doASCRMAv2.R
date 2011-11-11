##########################################################################
# Data set:
# ...
##########################################################################
library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "HapMap270,6.0,CEU,testSet";
chipType <- "GenomeWideSNP_6,Full";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doCRMAv2() on an AffymetrixCelSet object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doASCRMAv2() 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(dataSet, chipType=chipType, verbose=log);
print(res);
