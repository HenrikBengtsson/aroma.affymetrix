library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSet <- "GSE8605"
chipType <- "Mapping10K_Xba142";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
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
# RLE & NUSE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
qa <- QualityAssessmentModel(plm);
print(qa);

subplots(4, ncol=2, byrow=FALSE);
plotRle(qa);
plotNuse(qa);
# Reordered along x-axis
plotRle(qa, arrays=rev(seq(csR)));
plotNuse(qa, arrays=rev(seq(csR)));
