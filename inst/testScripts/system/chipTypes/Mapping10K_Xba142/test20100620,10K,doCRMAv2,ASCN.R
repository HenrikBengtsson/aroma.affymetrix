##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doCRMAv2()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- doCRMAv2(dataSet, chipType=chipType, combineAlleles=FALSE, verbose=log);
print(dsList);
