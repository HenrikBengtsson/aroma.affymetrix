###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the RMA 
# chip-effect estimates as estimated by affyPLM.
#
# Author: Mark Robinson and Henrik Bengtsson
# Created: 2007-06-20
# Last modified: 2008-07-30
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

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

csR <- AffymetrixCelSet$byName("Affymetrix-HeartBrain", 
                              chipType="HG-U133_Plus_2");
print(csR);

ab <- extractAffyBatch(csR, verbose=verbose);
print(ab);


###########################################################################
# HISTORY:
# 2009-10-16
# o Created.
###########################################################################
