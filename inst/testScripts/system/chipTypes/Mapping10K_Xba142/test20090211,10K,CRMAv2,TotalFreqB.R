##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE8605"
tags <- "ACC,-XY,v2,BPN,-XY,RMA,FLN,-XY";
chipType <- "Mapping10K_Xba142";
cesN <- SnpChipEffectSet$byName(dataSet, tags=tags, chipType=chipType, mergeStrands=TRUE);
print(cesN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Export to (total,fracB) ASB files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
theta <- getAromaUnitTotalCnBinarySet(cesN, verbose=log);
print(theta);

fracB <- getAromaUnitFracBCnBinarySet(cesN, verbose=log);
print(fracB);

units <- 100 + 1:10
print(extractMatrix(theta, units=units))
print(extractMatrix(fracB, units=units))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate average theta across arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thetaR <- calculateAverageColumnAcrossFiles(theta, verbose=log);
str(thetaR);

