##########################################################################
# Allele-specific CRMAv2 and Paired PSCBS
#
# Author: Henrik Bengtsson
# Created on: 2012-06-15
# Last updated: 2012-06-15
#
# DATA SET:
# GEO data set 'GSE12702'. Affymetrix CEL files are available from:
# 
#   http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12702
#
# Place them in rawData/GSE12702/Mapping250K_Nsp/*.CEL.  In total there
# are 20 tumor-normal pairs (40 CEL files).
##########################################################################
library("aroma.affymetrix");
library("calmate");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
extractSignals <- function(dsList, sampleName, reference=c("none", "median"), refIdxs=NULL, ..., verbose=FALSE) {
  reference <- match.arg(reference);
  idx <- indexOf(dsList$total, sampleName);
  dfT <- getFile(dsList$total, idx);
  dfB <- getFile(dsList$fracB, idx);
  tcn <- extractRawCopyNumbers(dfT, logBase=NULL, ..., verbose=verbose);
  baf <- extractRawAlleleBFractions(dfB, ..., verbose=verbose);
  if (reference == "median") {
    if (!is.null(refIdxs)) {
      dsR <- extract(dsList$total, refIdxs);
    } else {
      dsR <- dsList$total;
    }
    dfTR <- getAverageFile(dsR, verbose=verbose);
    tcnR <- extractRawCopyNumbers(dfTR, logBase=NULL, ..., verbose=verbose);
    tcn <- divideBy(tcn, tcnR);
    setSignals(tcn, 2*getSignals(tcn));
  }
  list(tcn=tcn, baf=baf);
} # extractSignals()


log <- Arguments$getVerbose(-8, timestamp=TRUE);

# Graphics options
devOptions("png", width=840);

dataSet <- "GSE12702";
chipType <- "Mapping250K_Nsp";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doASCRMAv2() on a single tumor-normal pair
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load all samples
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);

# Allele-specific CRMAv2
dsList <- doASCRMAv2(csR, verbose=log);
print(dsList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CalMaTe
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cmt <- CalMaTeCalibration(dsList);
print(cmt);

dsCList <- process(cmt, verbose=verbose);
print(dsCList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
