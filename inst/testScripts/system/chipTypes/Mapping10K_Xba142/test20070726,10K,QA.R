library("aroma.affymetrix")
log <- Verbose(threshold=-4, timestamp=TRUE);


dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, verbose=log);
keep <- 1:6;
csR <- extract(csR, keep);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fitting log-additive model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaPlm(csR);
print(plm);
fit(plm, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Some basic quality scores
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- getChipEffectSet(plm);

# Boxplots of log2(theta), RLE, and NUSE
layout(matrix(1:4, ncol=2, byrow=TRUE));
plotBoxplot(ces, type="theta", transform=log2);
plotBoxplot(ces, type="RLE", arrays=c(2,4:6));
plotBoxplot(ces, type="NUSE");
devDone();


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Advanced usage
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculating statistics without plotting
theta <- boxplotStats(ces, type="theta", transform=log2);
nuse <- boxplotStats(ces, type="NUSE");
rle <- boxplotStats(ces, type="RLE");

# Subset of arrays (avoids calculating stats for all arrays)
rleB <- boxplotStats(ces, type="RLE", arrays=c(2,4:6));
for (name in names(rleB)) {
  stopifnot(identical(rleB[name], rle[name]));
}

# Plotting the above statistics
layout(matrix(1:4, ncol=2, byrow=TRUE));
plotBoxplotStats(theta, main="theta");
plotBoxplotStats(rle[c(2,4:6)], main="RLE");
plotBoxplotStats(nuse, main="NUSE");
devDone();


# Calculates unit-specific RLE and NUSE scores
units <- 1000+1:5000;
theta <- extractMatrix(ces, units=units);
rle <- extractMatrix(ces, units=units, field="RLE");
nuse <- extractMatrix(ces, units=units, field="NUSE");

# Plotting the above statistics
layout(matrix(1:4, ncol=2, byrow=TRUE));
plotDensity(log2(theta), main="theta");
plotDensity(rle[,c(2,4:6)], main="RLE");
plotDensity(nuse, main="NUSE");
devDone();

# ...same, but basic unit annotation data added
units <- 1000+1:500;
theta <- extractDataFrame(ces, units=units, addNames=TRUE);
rle <- extractDataFrame(ces, units=units, field="RLE", addNames=TRUE);
nuse <- extractDataFrame(ces, units=units, field="NUSE", addNames=TRUE);
