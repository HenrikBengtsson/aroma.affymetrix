library("aroma.affymetrix");
library("ACNE");
log <- Arguments$getVerbose(-4, timestamp=TRUE);



dataSetName <- "HapMap270,100K,CEU,testSet";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

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
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csCList;
cesCnList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  plm <- NmfSnpPlm(cs, mergeStrands=TRUE);
  print(plm);
  fit(plm, verbose=log);
  ces <- getChipEffectSet(plm);
  print(ces);
  stopifnot(identical(getNames(ces), getNames(cs)));
  cesCnList[[chipType]] <- ces;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extraction test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thetaList <- list();
cesList <- cesCnList;
for (chipType in names(cesList)) {
  ces <- cesList[[chipType]];
  theta <- extractMatrix(ces, verbose=log);
  print(summary(theta));
  thetaList[[chipType]] <- theta;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesNList <- list();
for (chipType in names(csList)) {
  ces <- cesCnList[[chipType]];
  fln <- FragmentLengthNormalization(ces);
  print(fln);
  cesN <- process(fln, verbose=log);
  print(cesN);
  stopifnot(identical(getNames(cesN), getNames(ces)));
  cesNList[[chipType]] <- cesN;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Export total and fracB signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsTBList <- lapply(cesNList, FUN=function(cesN) {
  exportTotalAndFracB(cesN, verbose=verbose);
});

# Keep only the total CN data sets
dsTList <-lapply(dsTBList, FUN=function(dsList) dsList$total);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation model test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
seg <- CbsModel(dsTList);
print(seg);

fit(seg, arrays=1, chromosomes=19, verbose=log);

# Tests the case where one of the set does not have observations.
fit(seg, arrays=nbrOfArrays(seg), chromosomes=19, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ChromosomeExplorer test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ce <- ChromosomeExplorer(seg);
print(ce);
process(ce, arrays=1:2, chromosomes=c(19,23), verbose=log);

