library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);



dataSetName <- "HapMap270,100K,CEU,testSet";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993", 
                 "NA06994", "NA07000", "NA07019");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define a CEL set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipTypes[1], verbose=log);
keep <- 1:3;
cs <- extract(cs, keep);
sampleNames <- sampleNames[keep];
print(cs);
stopifnot(identical(getNames(cs), sampleNames));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set a custom CDF (here just use the other CDF of the chip in the set)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipTypes[2]);
setCdf(cs, cdf);
print(cs);
stopifnot(identical(getCdf(cs), cdf));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RMA background correction with custom CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bg <- RmaBackgroundCorrection(cs);
print(bg);
csBg <- process(bg, verbose=log);
print(csBg);
stopifnot(identical(getCdf(csBg), cdf));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Quantile normalization with custom CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
qn <- QuantileNormalization(cs);
print(qn);
csN <- process(qn, verbose=log);
print(csN);
stopifnot(identical(getCdf(csN), cdf));

