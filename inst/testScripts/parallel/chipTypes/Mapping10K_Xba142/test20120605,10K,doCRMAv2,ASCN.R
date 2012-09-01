##########################################################################
# Test script for running doCRMAv2() in parallel.  Note that this
# is still dangerous to do, because there may be rase conditions for:
#
#  - for creating/writing/updating the monocell CDF
#  - for writing cache files
#
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";

dsR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(dsR);

# Before running in parallel, assert that the monocell CDF has been created
cdf <- getCdf(dsR);
cdfM <- getMonoCell(cdf, verbose=verbose);

library("parallel");
cl <- makeCluster(2L);

clusterExport(cl, "dsR");
clusterExport(cl, "verbose");

res <- parLapply(cl, X=seq(dsR), fun=function(ii) {
  library("aroma.affymetrix");
  verbose && enter(verbose, sprintf("Array #%d", ii));
  dsListII <- doASCRMAv2(dsR, arrays=ii);
  verbose && exit(verbose);
  dsListII;
});

# Load result
dsNList <- doASCRMAv2(dsR);
print(dsNList);

