##########################################################################
# Allele-specific CRMAv2 and Paired PSCBS
#
# Author: Henrik Bengtsson
# Created on: 2011-11-11
# Last updated: 2012-09-02
##########################################################################
library("aroma.affymetrix");
library("PSCBS");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# Graphics options
devOptions("png", width=840);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE12702";
chipType <- "Mapping250K_Nsp";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract a single tumor-normal pair
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract tumor-normal pair
pair <- c(T="GSM318736", N="GSM318737");
csR <- extract(csR, indexOf(csR, pair));
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- doASCRMAv2(csR, verbose=verbose);
print(dsList);


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
fit <- segmentByPairedPSCBS(dfT, verbose=verbose);
print(fit);

# Plot segmentation results (saved to figures/)
pairName <- paste(pair, collapse="vs");
chrTag <- sprintf("Chr%s", seqToHumanReadable(getChromosomes(fit)));
toPNG(pairName, tags=c(chrTag, "PairedPSCBS"), aspectRatio=0.6, {
  plotTracks(fit);
});
