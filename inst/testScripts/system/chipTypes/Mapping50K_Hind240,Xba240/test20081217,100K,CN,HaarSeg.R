library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);



dataSetName <- "HapMap270,100K,CEU,testSet";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csRList <- list();
for (chipType in chipTypes) {
  cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
  print(cs);
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
  csCList[[chipType]] <- csC;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csCList;
cesCnList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  plm <- RmaCnPlm(cs, mergeStrands=TRUE, combineAlleles=TRUE, shift=300);
  print(plm);
  fit(plm, ram=1/2, verbose=log);
  ces <- getChipEffectSet(plm);
  print(ces);
  cesCnList[[chipType]] <- ces;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesCnList <- cesCnList;
cesNList <- list();
for (chipType in names(csList)) {
  ces <- cesCnList[[chipType]];
  fln <- FragmentLengthNormalization(ces);
  print(fln);
  cesN <- process(fln, verbose=log);
  print(cesN);
  cesNList[[chipType]] <- cesN;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Glad model test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
seg <- HaarSegModel(cesNList);
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

