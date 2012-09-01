library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSet <- "GSE9890";
chipType <- "HG-U133_Plus_2";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);

# GCRMA background correction requires an ACS file
cdf <- getCdf(csR);
acs <- getAromaCellSequenceFile(cdf);
print(acs);

bg <- GcRmaBackgroundCorrection(csR);
print(bg);
csB <- process(bg, verbose=verbose);
print(csB);

