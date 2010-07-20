library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "HapMap270"
tags <- "6.0,CEU,testSet,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY";
chipType <- "GenomeWideSNP_6,Full";
cesN <- SnpChipEffectSet$byName(dataSet, tags=tags, chipType=chipType, mergeStrands=TRUE);
print(cesN);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GC-content effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
unf <- getUnitNamesFile(cesN);
ugp <- getAromaUgpFile(unf);
print(ugp);
units <- getUnitsOnChromosomes(ugp, 1:22);
units <- units[seq(from=1, to=length(units), length.out=500e3)];

gcn <- GcContentNormalization2(cesN, target="zero");
devSet("GCN,before");
par(mar=c(4,3,1,1)+0.1, mgp=c(2.0,0.7,0));
plotCovariateEffects(gcn, units=units, ylim=c(0,16), verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GC-content normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
gcn <- GcContentNormalization2(cesN, target="zero");
print(gcn);
cesN2 <- process(gcn, verbose=log);
print(cesN2);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GC-content effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("GCN,after");
gcn2 <- GcContentNormalization2(cesN2, target="zero");
plotCovariateEffects(gcn2, units=units, ylim=c(0,16), verbose=log);

