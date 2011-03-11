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
# Quantile normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
qn <- QuantileNormalization(csR);
print(qn);
csN <- process(qn, verbose=log);
print(csN);
stopifnot(identical(getNames(csN), getNames(csR)));


