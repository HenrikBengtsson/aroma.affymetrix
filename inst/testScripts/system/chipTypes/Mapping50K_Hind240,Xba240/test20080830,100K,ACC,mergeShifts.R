library("aroma.affymetrix");
log <- Arguments$getVerbose(-4, timestamp=TRUE);



dataSetName <- "HapMap270,100K,CEU,testSet";
chipType <- "Mapping50K_Hind240";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR, mergeShifts=FALSE);
print(acc);
csC <- process(acc, verbose=log);
print(csC);
