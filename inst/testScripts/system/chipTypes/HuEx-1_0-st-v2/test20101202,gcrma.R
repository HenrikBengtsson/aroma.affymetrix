library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-3, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "Affymetrix-HeartBrain";
chipType <- "HuEx-1_0-st-v2";

cdf <- AffymetrixCdfFile$byChipType(chipType, tags="coreR3,A20071112,EP");
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="fullR3,A20071112,EP");
print(cdf);

# Setup CEL set using the core CDF.
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Background correction and normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdfS <- AffymetrixCdfFile$byChipType(getChipType(cdf, fullname=FALSE));
setCdf(csR, cdfS);
bc <- GcRmaBackgroundCorrection(csR, type="affinities");
print(bc);
csB <- process(bc, verbose=verbose);
print(csB);
setCdf(csB, cdf);
