##########################################################################
# Allele-specific CRMAv2
##########################################################################
library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE20584";
chipType <- "GenomeWideSNP_6,Full";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsNList <- doASCRMAv2(csR, verbose=verbose);
print(dsNList);

dsN <- exportAromaUnitPscnBinarySet(dsNList);
print(dsN);


db <- TabularTextFile("GSE20584,samples.txt");
setColumnNameTranslator(db, function(names, ...) {
  names <- gsub("id", "fixed", names);
  names <- gsub("fullname", "replacement", names);
  names;
});
df <- readDataFrame(db, colClassPatterns=c("*"="character"));
setFullNamesTranslator(dsN, df);