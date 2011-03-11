###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the estimates
# of the MAT (Model-based Analysis of Tiling arrays) algorithm.
#
# Author: Mark Robinson (and Henrik Bengtsson)
# Created: 2008-12-09
# Last modified: 2009-06-28
#
# Data set:
#  rawData/
#    MATtest/
#      Hs_PromPR_v02/
#        Prec1_MeDNA_IP1.CEL
#        Prec1_MeDNA_IP2.CEL
#        Prec1_MeDNA_Input1.CEL
#       TESTIP.bar.txt (estimates by MAT)
###########################################################################
library("aroma.affymetrix");
library("limma");  # makeContrast()
log <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup the tiling array data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02");
print(cdf);

csR <- AffymetrixCelSet$byName("MATtest", cdf=cdf);
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize the data using the MAT model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mn <- MatNormalization(csR);
csM <- process(mn, verbose=more(log, 3));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert data set such that it maps to the "unique" CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csU <- convertToUnique(csM, verbose=log);
print(csU);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Run 2 variations of MatSmoothing that will be compared to external estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sampleNames <- getNames(csU);

design1 <- makeContrasts(Prec1_MeDNA_IP1-Prec1_MeDNA_Input1, levels=sampleNames);
colnames(design1) <- "Prec1_IP1_minus_Input";
print(design1);

ms1 <- MatSmoothing(csU, design=design1, probeWindow=600, tags="singleIP");
csMS1 <- process(ms1, units=NULL, verbose=log);
print(csMS1);
stopifnot(nbrOfFiles(csMS1) == ncol(design1));

writeSgr(csMS1, verbose=log);


design2 <- makeContrasts(Prec1_MeDNA_IP1 + Prec1_MeDNA_IP2-Prec1_MeDNA_Input1, levels=sampleNames);
colnames(design2) <- "Prec1_IPs_minus_Input";

ms2 <- MatSmoothing(csU, design=design2, probeWindow=800, tags="multipleIP")
csMS2 <- process(ms2, units=NULL, verbose=log)
print(csMS2);
stopifnot(nbrOfFiles(csMS2) == ncol(design2));

writeSgr(csMS2, verbose=log);
