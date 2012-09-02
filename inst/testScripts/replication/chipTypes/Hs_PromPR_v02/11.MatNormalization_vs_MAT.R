###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the estimates
# of the MAT (Model-based Analysis of Tiling arrays) algorithm.
#
# Author: Mark Robinson and Henrik Bengtsson
# Created: 2008-12-09
# Last modified: 2012-09-02
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
#    GSE24546,testset/
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
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE24546";
tags <- "testset";
chipType <- "Hs_PromPR_v02";
sampleNamesMap <- c(
  GSM605951="Prec1_MeDNA_IP1",
  GSM605952="Prec1_MeDNA_IP2",
  GSM605953="Prec1_MeDNA_Input1"
);

cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

csR <- AffymetrixCelSet$byName(dataSet, tags=tags, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize the data using the MAT model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mn <- MatNormalization(csR);
csM <- process(mn, verbose=more(verbose, 3));
print(csM);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert data set such that it maps to the "unique" CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csU <- convertToUnique(csM, verbose=verbose);
print(csU);

# Rename
setFullNamesTranslator(csU, function(names, ...) { sampleNamesMap[names] });


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare to precomputed estimates from external MAT software
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare common units with prefix "chr1F".
cdfU <- getCdf(csU);
units <- indexOf(cdfU, "^chr1F");
cells <- getCellIndices(cdfU, units=units, stratifyBy="pm", 
                        unlist=TRUE, useNames=FALSE); 

# Get the chromosomal positions of these cells
acp <- AromaCellPositionFile$byChipType(getChipType(cdfU));
pos <- acp[cells,2,drop=TRUE]; 

# Order cells by chromsomal position
o <- order(pos);
pos <- pos[o];
cells <- cells[o];

# Extract the corresponding signals of the first array for each test
cf <- getFile(csU, 1);
print(cf);
sampleName <- getName(cf);

yN <- extractMatrix(cf, cells=cells, drop=TRUE, verbose=verbose);
str(yN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare to the MAT software
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAT normalized signals
normFile <- "Prec1_MeDNA_600_IP1-Input.tsv";
tsv <- TabularTextFile(filename=normFile, path=getPath(csR));
print(tsv);

# Read signals
colClassPatterns <- c("Chr"="character", "Pos"="integer", "numeric");
names(colClassPatterns)[3] <- sampleName;
data <- readDataFrame(tsv, colClassPatterns=colClassPatterns, nrow=435000);
data <- subset(data, (Chr == "chr1" & Pos %in% pos));

# Order as (pos,yN)
o <- match(data$Pos, pos);
# Sanity check
stopifnot(all(is.finite(o)));
data <- data[o,,drop=FALSE];

# Extract signals
yNB <- data[,3L];

# Sanity check
stopifnot(length(yNB) == length(yN));

# Compare on the log2 scale
yN <- log2(yN);

toPNG(getFullName(csU), tags=c(sampleName, "MatNormalization_vs_MAT"), {
  xlab <- expression(log[2](y[MAT]));
  ylab <- expression(log[2](y));
  plot(yNB, yN, pch=".", xlab=xlab, ylab=ylab);
  abline(a=0, b=1, col="red", lwd=2, lty=3);
  stext(side=3, pos=0, sampleName);
  stext(side=3, pos=1, "Chr01");
  stext(side=4, pos=1, sprintf("n=%d", length(yN)));
  stext(side=4, pos=0, getFullName(csU));
});

avgDiff <- median((yN-yNB)^2);
cat("Average difference: ", avgDiff, "\n", sep="");
stopifnot(avgDiff < 0.001);
