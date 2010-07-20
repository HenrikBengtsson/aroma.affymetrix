library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "HapMap270"
tags <- "6.0,CEU,testSet,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY";
chipType <- "GenomeWideSNP_6,Full";
cesN <- SnpChipEffectSet$byName(dataSet, tags=tags, chipType=chipType, mergeStrands=TRUE);
print(cesN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Export to (total,fracB) ASB files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- exportTotalAndFracB(cesN, verbose=log);

theta <- getAromaUnitTotalCnBinarySet(cesN, verbose=log);
print(theta);

fracB <- getAromaUnitFracBCnBinarySet(cesN, verbose=log);
print(fracB);

units <- 5000 + 1:10
print(extractMatrix(theta, units=units))
print(extractMatrix(fracB, units=units))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate average theta across arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thetaR <- calculateAverageColumnAcrossFiles(theta, verbose=log);
str(thetaR);
