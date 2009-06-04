library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);



dataSetName <- "Jeremy_2007-10k";
chipType <- "Mapping10K_Xba142";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
keep <- 1:6;
cs <- extract(cs, keep);
print(cs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(cs);
print(acc);
csC <- process(acc, verbose=log);
print(csC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=TRUE, shift=300);
print(plm);

fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces);
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Glad model test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
seg <- HaarSegModel(cesN);
print(seg);
fit(seg, arrays=1, chromosomes=19, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ChromosomeExplorer test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ce <- ChromosomeExplorer(seg, zooms=2^(0:3));
print(ce);
process(ce, chromosomes=c(19,23), verbose=log);
## process(ce, verbose=log);
