###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the RMA 
# background correction methods of affy and limma and that they
# generate the same results.
#
# Author: Henrik Bengtsson
# Created: 2010-10-29
# Last modified: 2010-10-29
#
# Data set:
#  rawData/
#   Affymetrix-HeartBrain/
#    HG-U133_Plus_2/
#     u1332plus_ivt_cerebellum_A.CEL [13555904 bytes]
#     u1332plus_ivt_cerebellum_B.CEL [13550687 bytes]
#     u1332plus_ivt_cerebellum_C.CEL [13551860 bytes]
#     u1332plus_ivt_heart_A.CEL      [13554731 bytes]
#     u1332plus_ivt_heart_B.CEL      [13553255 bytes]
#     u1332plus_ivt_heart_C.CEL      [13551203 bytes]
#  Source: Affymetrix Tissue samples, 2007.  http://www.affymetrix.com/
#  support/technical/sample_data/hugene_1_0_array_data.affx
###########################################################################

library("aroma.affymetrix");
library("affy");
library("limma");
log <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSet <- "Affymetrix-HeartBrain";
chipType <- "HG-U133_Plus_2";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# Using affy::bg.adjust()
bg <- RmaBackgroundCorrection(csR);
print(bg);
csB1 <- process(bg, verbose=log);
print(csB1);


# Using limma::backgroundCorrect()
bg <- NormExpBackgroundCorrection(csR);
print(bg);
csB2 <- process(bg, verbose=log);
print(csB2);

Y1 <- extractMatrix(csB1);
Y2 <- extractMatrix(csB2);
stopifnot(identical(Y1, Y2));



###########################################################################
# HISTORY:
# 2010-10-29 [HB]
# o Created from a private email 'NormExpBackgroundCorrection: Beta 
#   version available' on 2009-04-16.
###########################################################################
