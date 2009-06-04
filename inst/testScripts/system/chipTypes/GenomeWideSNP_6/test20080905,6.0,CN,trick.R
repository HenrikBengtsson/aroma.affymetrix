library("aroma.affymetrix");

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

if (length(findUnitsTodo(plm)) > 0) {
   # Fit CN probes quickly (~5-10s/array + some overhead)
  units <- fitCnProbes(plm, verbose=log);
  str(units);
  # int [1:945826] 935590 935591 935592 935593 935594 935595 ...

  # Fit remaining units, i.e. SNPs (~5-10min/array)
  units <- fit(plm, verbose=log);
  str(units);
}

ces <- getChipEffectSet(plm);
print(ces);

theta <- extractMatrix(ces, units=1000:1002);

fln <- FragmentLengthNormalization(ces);
cesN <- process(fln, verbose=log);
print(cesN);

cnr <- CbsModel(cesN);
print(cnr);

ce <- ChromosomeExplorer(cnr);
print(ce);
process(ce, arrays=1, chromosomes=c(19,23), verbose=log);
