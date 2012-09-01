##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE8605"
tags <- "ACC,-XY,BPN,-XY,RMA,+300,FLN,-XY";
chipType <- "Mapping10K_Xba142";
cesN <- SnpChipEffectSet$byName(dataSet, tags=tags, chipType=chipType, mergeStrands=TRUE);
print(cesN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Export to (total,fracB) ASB files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
theta <- getAromaUnitTotalCnBinarySet(cesN, verbose=verbose);
print(theta);

fracB <- getAromaUnitFracBCnBinarySet(cesN, verbose=verbose);
print(fracB);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate average theta across arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thetaR <- calculateAverageColumnAcrossFiles(theta, verbose=verbose);
str(thetaR);

