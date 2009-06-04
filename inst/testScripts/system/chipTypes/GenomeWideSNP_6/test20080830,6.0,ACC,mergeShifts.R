library("aroma.affymetrix");

log <- Arguments$getVerbose(-50, timestamp=TRUE);

dataSetName <- "HapMap270,6.0,CEU,testSet";
chipType <- "GenomeWideSNP_6";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up CEL set and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full");
print(cdf);

csR <- AffymetrixCelSet$byName(dataSetName, cdf=cdf, verbose=log);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic-crosstalk calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR, mergeShifts=FALSE);
print(acc);

csC <- process(acc, verbose=log);
print(csC);
