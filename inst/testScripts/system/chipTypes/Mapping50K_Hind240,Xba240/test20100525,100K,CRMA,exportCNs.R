library("aroma.affymetrix")

verbose <- Verbose(threshold=-4, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doCBS() with explicit data set tuple
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap270,100K,CEU,testSet";
tags <- "ACC,-XY,RMA,+300,A+B,FLN,-XY";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

nbrOfSets <- length(chipTypes);
for (kk in seq(length=nbrOfSets)) {
  chipType <- chipTypes[kk];
  ds <- CnChipEffectSet$byName(dataSet, tags=tags, chipType=chipType,
                              mergeStrands=TRUE, combineAlleles=TRUE);

  dsOut <- exportTotalAndFracB(ds, verbose=verbose);
  verbose && print(verbose, dsOut);
}


