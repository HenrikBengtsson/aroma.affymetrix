library("aroma.affymetrix");
log <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup the tiling array data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02", tags="Harvard,ROIs");
print(cdf);

csR <- AffymetrixCelSet$byName("MNtest", cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize the data using the MAT model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mn <- MatNormalization(csR, numChunks=20);
csM <- process(mn, verbose=more(log, 30));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert data set such that it maps to the "unique" CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csU <- convertToUnique(csM, verbose=log);
print(csU);

dm <- matrix(1, nrow=nbrOfFiles(csU), ncol=3);
rownames(dm) <- getNames(csU);
colnames(dm) <- sprintf("MATSmoothingTest-%d", 1:ncol(dm));
ms <- MatSmoothing(csU, design=dm, probeWindow=300, nProbes=10);
print(ms);

dsMS <- process(ms, units=1:100, verbose=verbose);
print(dsMS);
