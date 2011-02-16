library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);


## dataSet <- "MNtest";
## chipType <- "Hs_PromPR_v02,Harvard,ROIs";

dataSet <- "E-MEXP-1481";
chipType <- "Hs_PromPR_v02"; 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup the tiling array data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize the data using the MAT model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mn <- MatNormalization(csR);
csM <- process(mn, verbose=more(verbose, 30));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert data set such that it maps to the "unique" CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csU <- convertToUnique(csM, verbose=verbose);
print(csU);

dm <- matrix(1, nrow=nbrOfFiles(csU), ncol=3);
rownames(dm) <- getNames(csU);
colnames(dm) <- sprintf("MATSmoothingTest-%d", 1:ncol(dm));
ms <- MatSmoothing(csU, design=dm, probeWindow=300, nProbes=10);
print(ms);

dsMS <- process(ms, units=1:100, verbose=verbose);
print(dsMS);
