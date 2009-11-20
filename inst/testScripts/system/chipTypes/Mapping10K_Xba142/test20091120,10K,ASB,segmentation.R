##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE8605"
tags <- "ACC,-XY,v2,BPN,-XY,RMA,FLN,-XY";
chipType <- "Mapping10K_Xba142";

ds <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
print(ds);

cbs <- CbsModel(ds);
print(cbs);

library("aroma.affymetrix");  # AD HOC fix to find UGP file. /HB 2009-11-20
fit(cbs, arrays=1:2, chromosomes=5:8, verbose=verbose);

ce <- ChromosomeExplorer(cbs);
print(ce);

process(ce, arrays=1:2, chromosomes=1:8, verbose=verbose);
