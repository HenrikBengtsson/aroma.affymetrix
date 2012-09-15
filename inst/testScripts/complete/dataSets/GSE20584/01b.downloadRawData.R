path <- system.file("testScripts", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");



##########################################################################
# Data set:
# GSE20584
#   GenomeWideSNP_6/
#    *.CEL [2]
#
# Overall design:
#  One lung tumor sample and an adjacent normal sample were assayed on
#  the Affymetrix SNP6.0 array.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE20584
##########################################################################
dataSet <- "GSE20584";
chipType <- "GenomeWideSNP_6";

verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType, 
                   chipTypeAliases=c("GenomeWideEx_6"="GenomeWideSNP_6"));
print(ds);


verbose && exit(verbose);
