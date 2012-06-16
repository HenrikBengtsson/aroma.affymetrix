##########################################################################
# Allele-specific CRMAv2 and Paired PSCBS
#
# Author: Henrik Bengtsson
# Created on: 2011-11-11
# Last updated: 2011-11-11
#
# DATA SET:
# GEO data set 'GSE12702'. Affymetrix CEL files are available from:
# 
#   http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12702
#
# Place them in rawData/GSE12702/Mapping250K_Nsp/*.CEL.  In total there
# are 20 tumor-normal pairs (40 CEL files).  However, for this example
# we only need two of the CEL files, namely GSM318736 (tumor) and
# GSM318737 (matched normal).
##########################################################################
library("aroma.affymetrix");
library("PSCBS");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

# Graphics options
devOptions("png", width=840);

dataSet <- "GSE12702";
chipType <- "Mapping250K_Nsp";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doASCRMAv2() on a single tumor-normal pair
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load all samples
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);

# Extract tumor-normal pair
pair <- c(T="GSM318736", N="GSM318737");
csR <- extract(csR, indexOf(csR, pair));
print(csR);

# Allele-specific CRMAv2 on pair
res <- doASCRMAv2(csR, verbose=log);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get (chromosome, position) annotation data
ugp <- getAromaUgpFile(res$total);
chromosome <- ugp[,1,drop=TRUE];
x <- ugp[,2,drop=TRUE];

# Extract (total,beta) estimates for the tumor-normal pair
data <- extractPSCNArray(res$total);
dimnames(data)[[3]] <- names(pair);
str(data);

CT <- 2 * (data[,"total","T"] / data[,"total","N"]);
betaT <- data[,"fracB","T"];
betaN <- data[,"fracB","N"];

# Setup data structure for Paired PSCBS
df <- data.frame(chromosome=chromosome, x=x, CT=CT, betaT=betaT, betaN=betaN);

# Segment Chr 8.
dfT <- subset(df, chromosome == 8);

# Paired PSCBS
fit <- segmentByPairedPSCBS(dfT, verbose=log);
print(fit);

# Plot segmentation results (saved to figures/)
pairName <- paste(pair, collapse="vs");
chrTag <- sprintf("Chr%s", seqToHumanReadable(getChromosomes(fit)));
toPNG(pairName, tags=c(chrTag, "PairedPSCBS"), aspectRatio=0.6, {
  plotTracks(fit);
});
