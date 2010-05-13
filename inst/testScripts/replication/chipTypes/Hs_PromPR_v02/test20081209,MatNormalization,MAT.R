###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the estimates
# of the MAT (Model-based Analysis of Tiling arrays) algorithm.
#
# Author: Mark Robinson (and Henrik Bengtsson)
# Created: 2008-12-09
# Last modified: 2009-09-07
#
# Data set:
#  annotationData/
#    chipTypes/
#      Hs_PromPR_v02/
#        Hs_PromPR_v02.cdf
#        Hs_PromPR_v02.acs
#        Hs_PromPR_v02.acm
#        Hs_PromPR_v02,unique.acp
#  rawData/
#    MATtest/
#      Hs_PromPR_v02/
#        Prec1_MeDNA_IP1.CEL
#        Prec1_MeDNA_IP2.CEL
#        Prec1_MeDNA_Input1.CEL
#        # Validation files
#        Prec1_MeDNA_600_IP1-Input.tsv
#        Prec1_MeDNA_600_IP1-Input.bar.txt
#        Prec1_MeDNA_800_IPs-Input.bar.txt
###########################################################################
library("aroma.affymetrix");
library("limma");  # makeContrast()
log <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert that required annotation data files exist
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02");
print(cdf);

acs <- AromaCellSequenceFile$byChipType(getChipType(cdf));
print(acs);

acm <- AromaCellMatchScoreFile$byChipType(getChipType(cdf));
print(acm);

cdfU <- getUniqueCdf(cdf, verbose=log);
print(cdfU);

acp <- AromaCellPositionFile$byChipType(getChipType(cdfU));
print(acp);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup the tiling array data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName("MATtest", cdf=cdf);
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize the data using the MAT model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mn <- MatNormalization(csR, numChunks=20);
csM <- process(mn, verbose=more(log, 3));
print(csM);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert data set such that it maps to the "unique" CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csU <- convertToUnique(csM, verbose=log);
print(csU);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Run 2 variations of MatSmoothing that will be compared to 
# external estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sampleNames <- getNames(csU);

design1 <- makeContrasts(Prec1_MeDNA_IP1-Prec1_MeDNA_Input1, levels=sampleNames);
colnames(design1) <- "Prec1_IP1_minus_Input";
print(design1);

ms1 <- MatSmoothing(csU, design=design1, probeWindow=600, tag="singleIP");
csMS1 <- process(ms1, units=NULL, verbose=log);
print(csMS1);
stopifnot(nbrOfFiles(csMS1) == ncol(design1));

design2 <- makeContrasts(Prec1_MeDNA_IP1 + Prec1_MeDNA_IP2-Prec1_MeDNA_Input1, levels=sampleNames);
colnames(design2) <- "Prec1_IPs_minus_Input";

ms2 <- MatSmoothing(csU, design=design2, probeWindow=800, tag="multipleIP")
csMS2 <- process(ms2, units=NULL,verbose=log)
print(csMS2);
stopifnot(nbrOfFiles(csMS2) == ncol(design2));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare to external estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare common units with prefix "chr1F".
cdf <- getCdf(csU);
units <- indexOf(cdf, "^chr1F");
cells <- getCellIndices(cdf, units=units, stratifyBy="pm", 
                        unlist=TRUE, useNames=FALSE); 

# Get the chromosomal positions of these cells
acp <- AromaCellPositionFile$byChipType(getChipType(cdf));
pos <- acp[cells,2,drop=TRUE]; 

# Order cells by chromsomal position
o <- order(pos);
pos <- pos[o];
cells <- cells[o];

# Extract the corresponding signals of the first array for each test
cf <- getFile(csU, 1);
normSampleName <- getName(cf);
yN <- extractMatrix(cf, cells=cells, drop=TRUE, verbose=log);

cf <- getFile(csMS1, 1);
y1 <- extractMatrix(cf, cells=cells, drop=TRUE, verbose=log);

cf <- getFile(csMS2, 1);
y2 <- extractMatrix(cf, cells=cells, drop=TRUE, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load external and compare normalized signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# file with normalized signals
normFile <- "Prec1_MeDNA_600_IP1-Input.tsv";

tsv <- TabularTextFile(filename=normFile, path=getPath(csR));
data <- readDataFrame(tsv, colClasses=colClasses("cinn"), nrow=435000);
data <- subset(data, (Chr == "chr1" & Pos %in% pos));
o <- order(data$Pos);
data <- data[o,,drop=FALSE];

col <- grep(normSampleName, getColumnNames(tsv));
yNB <- data[,col,drop=TRUE];

# Compare on the log2 scale
yN <- log2(yN);

plot(yN, yNB, pch=".");
abline(a=0, b=1, col="red", lwd=2, lty=3);  

stopifnot(length(yN) == length(yNB));
avgDiff <- median((yN-yNB)^2);
cat(avgDiff, "\n");
stopifnot(avgDiff < 0.001);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load external and compare smoothed normalized signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# files with native-MAT summarized signals
summaryFiles <- c("Prec1_MeDNA_600_IP1-Input.bar.txt", 
                  "Prec1_MeDNA_800_IPs-Input.bar.txt");

yRefList <- lapply(summaryFiles, FUN=function(filename) {
  columnNames <- c("chromosome", "position", "signal");
  bar <- TabularTextFile(filename, path=getPath(csR), columnNames=columnNames);
  data <- readDataFrame(bar, colClasses=colClasses("cin"), nrow=435000);

  # Extract the subset available in the aroma.affymetrix estimate
  data <- subset(data, chromosome == "chr1" & position %in% pos);

  # Order by position
  o <- order(data$position);
  data <- data[o,,drop=FALSE];
  
  # Extract the external estimates
  data$signal;
});

# Compare on the log2 scale
y1 <- log2(y1);
y2 <- log2(y2);

yList <- list(y1, y2);

for (ii in 1:length(yList)) {
  y <- yList[[ii]];
  yRef <- yRefList[[ii]];
  # Sanity check  
  stopifnot(length(y) == length(yRef));

  plot(y, yRef, pch=".");
  abline(a=0, b=1, col="red", lwd=2, lty=3);  

  # Our CDF scores a few probes that the external data doesn't, 
  # external data has these as 0
  keep <- which(yRef != 0);
  avgDiff <- median((y[keep]-yRef[keep])^2);
  cat(avgDiff, "\n");
  stopifnot(avgDiff < 0.06);
} # for (ii ...)

print(sessionInfo());
