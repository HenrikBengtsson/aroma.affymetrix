library("aroma.affymetrix");

ces <- doRMA("Affymetrix-HeartBrain", chipType="HG-U133_Plus_2", verbose=-5);
print(ces);

eset <- extractExpressionSet(ces, verbose=-5);
print(eset);
print(sampleNames(eset));

# Sanity checks
stopifnot(identical(sampleNames(eset), getNames(ces)));
stopifnot(identical(featureNames(eset), getUnitNames(getCdf(ces))));
