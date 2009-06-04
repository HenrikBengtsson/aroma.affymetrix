library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);



dataSetName <- "HapMap270,100K,CEU,testSet";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");
#chipTypes <- chipTypes[2];

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993", 
                 "NA06994", "NA07000", "NA07019");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csRList <- list();
for (chipType in chipTypes) {
  cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
  print(cs);
  stopifnot(identical(getNames(cs), sampleNames));
  csRList[[chipType]] <- cs;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csRList;
csCList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  acc <- AllelicCrosstalkCalibration(cs);
  print(acc);
  csC <- process(acc, verbose=log);
  print(csC);
  stopifnot(identical(getNames(csC), getNames(cs)));
  csCList[[chipType]] <- csC;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot allele pairs before and after calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (what in c("input", "output")) {
  plotAllelePairs(acc, array=1, what=what, verbose=log);
}
