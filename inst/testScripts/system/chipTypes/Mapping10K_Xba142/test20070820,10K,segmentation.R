library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, verbose=log);
keep <- 1:6;
csR <- extract(csR, keep);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR);
print(acc);
csC <- process(acc, verbose=log);
print(csC);
stopifnot(identical(getNames(csC), getNames(csR)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=TRUE, shift=300); 
print(plm);

fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);
stopifnot(identical(getNames(ces), getNames(csR)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces);
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);
stopifnot(identical(getNames(cesN), getNames(ces)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Try to setup different segmentation models
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
constList <- list(
  glad = GladModel,
  cbs  = CbsModel,
  haar = HaarSegModel
);
print(names(constList));

csmList <- lapply(constList, FUN=function(csm) {
  res <- NULL;
  tryCatch({
    res <- csm(cesN);
  }, error = function(ex) {
    print(ex);
  })
  res;
});

csmList <- csmList[!sapply(csmList, is.null)];
print(names(constList));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Try to segment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitList <- lapply(csmList, FUN=function(csm) {
  fit(csm, arrays=1, chromosomes=19, verbose=log);
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Try generate ChromosomeExplorer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
lapply(csmList, FUN=function(csm) {
  ce <- ChromosomeExplorer(csm);
  process(ce, arrays=1:2, chromosomes=c(1:2,23), verbose=log);
});

