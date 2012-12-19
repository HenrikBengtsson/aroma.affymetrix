library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE9890";
chipType <- "HG-U133_Plus_2";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);

ab <- extractAffyBatch(csR, verbose=verbose);
print(ab);

ces <- doRMA(csR, verbose=verbose);
print(ces);

eset <- extractExpressionSet(ces, verbose=verbose);
print(eset);
print(sampleNames(eset));

# Sanity checks
stopifnot(identical(sampleNames(eset), getNames(ces)));
stopifnot(identical(featureNames(eset), getUnitNames(getCdf(ces))));
