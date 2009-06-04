library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);



dataSetName <- "Jeremy_2007-10k";
chipType <- "Mapping10K_Xba142";

# Expected sample names
sampleNames <- c("0001-7", "0002-10", "0004-13", "0005-14", "0007-18", 
                      "0008-19", "0010-22", "2-DPrrr", "MH12", "MH18");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
keep <- 1:6;
csR <- extract(csR, keep);
sampleNames <- sampleNames[keep];
print(csR);
stopifnot(identical(getNames(csR), sampleNames));


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
