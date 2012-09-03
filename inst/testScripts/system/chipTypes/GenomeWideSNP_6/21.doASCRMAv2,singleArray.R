library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE13372,testset";
chipType <- "GenomeWideSNP_6,Full";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);

csR <- extract(csR, 1);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (a) AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (b) CRMAv2 - single array
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup a single array
csR1 <- extract(csR, subset[1]);

# Rename data set in order to not pick up existing results
setFullName(csR1, sprintf("%s,singleArray", getFullName(csR1)));

dsNList1 <- doASCRMAv2(csR1, verbose=verbose);
print(dsNList1);


# Sanity checks
dsNList0 <- lapply(dsNList, FUN=extract, 1);
for (key in names(dsNList1)) {
  dsN0 <- dsNList0[[key]];
  dsN1 <- dsNList1[[key]];
  stopifnot(getFullNames(dsN1) == getFullNames(dsN0));
  data0 <- extractMatrix(dsN0);
  data1 <- extractMatrix(dsN1);
  stopifnot(all.equal(data1, data0));
}
