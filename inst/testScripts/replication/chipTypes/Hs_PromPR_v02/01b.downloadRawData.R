path <- system.file("testScripts", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");



##########################################################################
# Data set:
# GSE24546,testset/
#   Hs_PromPR_v02/
#     *.CEL [3]
#
# Overall design:
#  Comparison of MeDIP/MBD for DNA methylation profiling, comparison of
#  whole genome amplification techniques, using tiling array for copy
#  number aberration detection and comparisons of tiling array data to
#  sequencing readouts.
#
# URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24546
##########################################################################
dataSet <- "GSE24546";
tags <- "testset";
chipType <- "Hs_PromPR_v02";
sampleNamesMap <- csampleNamesMap <- c(
  GSM605951="Prec1_MeDNA_IP1",
  GSM605952="Prec1_MeDNA_IP2",
  GSM605953="Prec1_MeDNA_Input1"
);

verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataFiles(dataSet, tags=tags, chipType=chipType, sampleNames=names(sampleNamesMap));
print(ds);
## AffymetrixCelSet:
## Name: GSE24546
## Tags:
## Path: rawData/GSE24546/Hs_PromPR_v02
## Platform: Affymetrix
## Chip type: Hs_PromPR_v02
## Number of arrays: 3
## Names: GSM605947, GSM605946, GSM605945 [3]
## Time period: 2007-10-30 12:51:27 -- 2007-10-30 13:59:56
## Total file size: 134.53MB
## RAM: 0.01MB

verbose && exit(verbose);
