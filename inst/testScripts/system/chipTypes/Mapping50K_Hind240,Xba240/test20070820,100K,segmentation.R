library("aroma.affymetrix")

log <- Verbose(threshold=-4, timestamp=TRUE);

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
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csCList;
cesCnList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  plm <- RmaCnPlm(cs, mergeStrands=TRUE, combineAlleles=TRUE, shift=300);
  print(plm);
  fit(plm, verbose=log);
  ces <- getChipEffectSet(plm);
  print(ces);
  stopifnot(identical(getNames(ces), getNames(cs)));
  cesCnList[[chipType]] <- ces;
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
# Set up list of CN segmenation models
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cnsList <- list(
  "glad" = GladModel(cesNList),
   "cbs" = CbsModel(cesNList)
);
print(cnsList);

# Fit a subset
lapply(cnsList, FUN=fit, arrays=1, chromosomes=19, verbose=log);

# There is a CN deletion on chr 2 @ 83.0Mb in NA06985.
lapply(cnsList, FUN=function(cns) {
  ce <- ChromosomeExplorer(cns);
  process(ce, arrays=1, chromosomes=c(2,19), verbose=log);
})

# There is a CN deletion on chr 22 @ 21Mb in NA06994.
lapply(cnsList, FUN=function(cns) {
  ce <- ChromosomeExplorer(cns);
  process(ce, arrays=4, chromosomes=22, verbose=log);
})

