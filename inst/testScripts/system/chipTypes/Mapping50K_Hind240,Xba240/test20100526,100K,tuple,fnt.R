library("aroma.affymetrix")

verbose <- Verbose(threshold=-4, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doCBS() with explicit data set tuple
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap270,100K,CEU,testSet";
tags <- "ACC,-XY,RMA,+300,A+B,FLN,-XY";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

nbrOfSets <- length(chipTypes);
dsList <- vector("list", nbrOfSets);
for (kk in seq(length=nbrOfSets)) {
  chipType <- chipTypes[kk];
  ds <- CnChipEffectSet$byName(dataSet, tags=tags, chipType=chipType,
                              mergeStrands=TRUE, combineAlleles=TRUE);

  fnts <- getAromaFullNameTranslatorSet(ds);
  print(fnts);
  appendFullNamesTranslator(ds, fnts);

  cat(verbose, "Default fullnames:");
  print(verbose, getFullNames(ds, translate=FALSE));
  cat(verbose, "Translated fullnames:");
  print(verbose, getFullNames(ds));

  dsList[[kk]] <- ds;
}
print(dsList);

dsTuple <- as.CopyNumberDataSetTuple(dsList);
print(dsTuple);
print(getNames(dsTuple));

