library("aroma.affymetrix")

log <- Verbose(threshold=-50, timestamp=TRUE);

dataSetName <- "HapMap270,6.0,CEU,testSet";
chipType <- "GenomeWideSNP_6";

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993", 
                 "NA06994", "NA07000", "NA07019");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up CEL set and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full");
print(cdf);

csR <- AffymetrixCelSet$byName(dataSetName, cdf=cdf, verbose=log);
print(csR);
stopifnot(identical(getNames(csR), sampleNames));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic-crosstalk calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR);
print(acc);

csC <- process(acc, verbose=log);
print(csC);
stopifnot(identical(getNames(csC), getNames(csR)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Averaging probe-level model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- AvgCnPlm(csC, mergeStrands=TRUE, combineAlleles=TRUE, shift=300);
print(plm);

ces <- getChipEffectSet(plm);
print(ces);

if (length(findUnitsTodo(plm)) > 0) {
   # Fit CN probes quickly (~5-10s/array + some overhead)
  units <- fitCnProbes(plm, verbose=log);
  str(units);
  # int [1:945826] 935590 935591 935592 935593 935594 935595 ...

  # Fit remaining units, i.e. SNPs (~5-10min/array)
  units <- fit(plm, verbose=log);
  str(units);
}


fln <- FragmentLengthNormalization(ces);
cesN <- process(fln, verbose=log);
print(cesN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation with specific reference set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Use the robust average of the first three arrays as a reference
cesR <- extract(cesN, 1:3);
ceR <- getAverageFile(cesR);
print(ceR);

sm <- CbsModel(cesN, ceR);
print(sm);
fit(sm, arrays=1, chromosomes=19, verbose=log);
