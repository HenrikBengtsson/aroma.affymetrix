library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);



dataSetName <- "HapMap270,100K,CEU,testSet";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");
chipTypes <- chipTypes[1];

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
# Quantile normalization tests #1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csNList <- list();
csList <- csRList;
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  qn <- QuantileNormalization(cs);
  print(qn);
  csN <- process(qn, verbose=log);
  print(csN);
  stopifnot(identical(getNames(csN), getNames(cs)));
  csNList[[chipType]] <- csN;
}
