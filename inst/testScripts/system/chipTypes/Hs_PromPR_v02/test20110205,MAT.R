library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);


dataSet <- "E-MEXP-1481";
chipType <- "Hs_PromPR_v02";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup the tiling array data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize the data using the MAT model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mn <- MatNormalization(csR);
csM <- process(mn, verbose=verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert data set such that it maps to the "unique" CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get the "unique" CDF, which is generated if missing
cdfU <- getUniqueCdf(cdf, verbose=verbose);
print(cdfU);

csU <- convertToUnique(csM, verbose=verbose);
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
y <- extractMatrix(cf, cells=cells, drop=TRUE, verbose=verbose);

par(mar=c(3,3,1,1)+0.1, mgp=c(1.8,0.5,0));
xlab <- "Physical position (Mb)";
ylab <- expression(log2(PM));
plot(pos/1e6, log2(y), pch=".", cex=2, xlab=xlab, ylab=ylab);
stext(side=3, pos=0, cex=0.7, getFullName(cf));
stext(side=3, pos=1, cex=0.7, sprintf("Chr%02d (n=%d)", chr, length(pos)));
stext(side=4, pos=1, cex=0.7, getChipType(cdf));
