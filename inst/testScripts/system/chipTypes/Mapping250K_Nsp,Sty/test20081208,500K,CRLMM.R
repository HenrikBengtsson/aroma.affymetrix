library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap270,500K,CEU,testSet";
chipType <- "Mapping250K_Nsp";

cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA according to aroma.affymetrix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- justSNPRMA(csR, normalizeToHapmap=TRUE, returnESet=FALSE, verbose=log);
print(ces);

crlmm <- CrlmmModel(ces);
print(crlmm);

units <- fit(crlmm, verbose=log);
str(units);

gcs <- getCallSet(crlmm);
print(gcs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot along genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chromosome <- 2;
gi <- getGenomeInformation(cdf);
units <- getUnitsOnChromosome(gi, chromosome=chromosome);
pos <- getPositions(gi, units=units)/1e6;

layout(matrix(seq(ces), ncol=2, byrow=TRUE));
par(mar=c(3.5,4,1.5,1), mgp=c(1.8,0.5,0), pch=".");
for (array in seq(ces)) {
  ce <- getFile(ces, array);
  gc <- getFile(gcs, array);
  data <- extractTotalAndFreqB(ce, units=units, drop=TRUE);
  calls <- extractGenotypes(gc, units=units, drop=TRUE);
  col <- c(AA=1, AB=2, BB=3)[calls];
  plot(pos, data[,"freqB"], col=col, cex=4, ylim=c(0,1));
  stext(side=3, pos=0, getName(ce));
  stext(side=3, pos=1, sprintf("Chr%02d", chromosome));
}
