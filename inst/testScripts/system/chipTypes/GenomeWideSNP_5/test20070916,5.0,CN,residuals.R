library("aroma.affymetrix")

log <- Arguments$getVerbose(-50, timestamp=TRUE);



dataSetName <- "HapMap270,5.0,CEU,testSet";
chipType <- "GenomeWideSNP_5";

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993",
                 "NA07019", "NA07022", "NA07056");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up CEL set and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full,r2");
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
# Probe-level model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaPlm(csC);
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Residuals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rs <- calculateResidualSet(plm, verbose=log);
print(rs);

ae <- ArrayExplorer(rs);
setColorMaps(ae, c("log2,log2neg,rainbow", "log2,log2pos,rainbow"));
print(ae);
process(ae, arrays=1:2, verbose=log);
