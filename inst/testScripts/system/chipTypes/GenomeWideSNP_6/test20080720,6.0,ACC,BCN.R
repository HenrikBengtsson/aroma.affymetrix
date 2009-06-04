library("aroma.affymetrix")

log <- Arguments$getVerbose(-50, timestamp=TRUE);

dataSetName <- "HapMap270,6.0,CEU,testSet";
chipType <- "GenomeWideSNP_6";

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993", 
                 "NA06994", "NA07000", "NA07019");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up CEL set and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full");
print(cdf);

csR <- AffymetrixCelSet$byName(dataSetName, cdf=cdf, verbose=log);
print(csR);
stopifnot(identical(getNames(csR), sampleNames));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic-crosstalk calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR);
print(acc);

csC <- process(acc, verbose=log);
print(csC);
stopifnot(identical(getNames(csC), getNames(csR)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Base-count normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bcn <- BaseCountNormalization(csC, shift=+300, bootstrap=TRUE);
print(bcn);

csN <- process(bcn, verbose=log);
print(csN);
stopifnot(identical(getNames(csN), getNames(csR)));
