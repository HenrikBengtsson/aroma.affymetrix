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
cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
keep <- 1:6;
cs <- extract(cs, keep);
sampleNames <- sampleNames[keep];
print(cs);
stopifnot(identical(getNames(cs), sampleNames));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(cs);
print(acc);
csC <- process(acc, verbose=log);
print(csC);
stopifnot(identical(getNames(csC), getNames(cs)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=TRUE, shift=300);
print(plm);

fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);
stopifnot(identical(getNames(ces), getNames(cs)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces);
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);
stopifnot(identical(getNames(cesN), getNames(ces)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation with specific reference set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Use the robust average of the first five arrays as a reference
cesR <- extract(cesN, 1:5);
ceR <- getAverageFile(cesR);
print(ceR);

sm <- CbsModel(cesN, ceR);
print(sm);
fit(sm, arrays=1, chromosomes=19, verbose=log);


