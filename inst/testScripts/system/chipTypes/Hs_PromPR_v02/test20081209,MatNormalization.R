library("aroma.affymetrix");
log <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup the tiling array data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02", tags="Harvard,ROIs");
print(cdf);

csR <- AffymetrixCelSet$byName("MNtest", cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize the data using the MAT model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mn <- MatNormalization(csR, numChunks=20);
csM <- process(mn, verbose=more(log, 30));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert data set such that it maps to the "unique" CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get the "unique" CDF, which is generated if missing
cdfU <- getUniqueCdf(cdf, verbose=more(log, 60));
print(cdfU);

csU <- convertToUnique(csM, verbose=log);
print(csU);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot some of the estimates on a single chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cs <- csU;
cdf <- getCdf(cs);
acp <- AromaCellPositionFile$byChipType(getChipType(cdf));
stopifnot(nbrOfCells(acp) == nbrOfCells(cdf));

# Identify cells on Chr21
chr <- 21;
cells <- whichVector(acp[,1,drop=TRUE] == chr);
str(cells);

# Get their chromosomal positions
pos <- acp[cells,2,drop=TRUE];

# Extract data for array #1
cf <- getFile(cs, 1);
y <- extractMatrix(cf, cells=cells, drop=TRUE, verbose=log);

par(mar=c(3,3,1,1)+0.1, mgp=c(1.8,0.5,0));
xlab <- "Physical position (Mb)";
ylab <- expression(log2(PM));
plot(pos/1e6, log2(y), pch=".", cex=2, xlab=xlab, ylab=ylab);
stext(side=3, pos=0, cex=0.7, getFullName(cf));
stext(side=3, pos=1, cex=0.7, sprintf("Chr%02d (n=%d)", chr, length(pos)));
stext(side=4, pos=1, cex=0.7, getChipType(cdf));
