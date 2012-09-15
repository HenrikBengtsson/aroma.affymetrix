path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# GSE34754/
#   Mapping250K_Nsp/
#    GSM854615.CEL, ..., GSM854626.CEL [12]
#
# Overall design:
#  Breast cancer cell lines HCC38 and HCC1143 and paired B 
#  lymphoblastoid cell lines HCC38-BL and HCC1143-BL were
#  purchased from ATCC. [...] To mimic tumor containing normal cells,
#  DNA from HCC38 and HCC1143 cells was mixed with DNA from autologous
#  B lymphoblastoid cells HCC38-BL and HCC1143-BL, respectively,
#  in ratios (w/w) 100:0, 80:20, 60:40, 40:60, and 20:80.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE34754
##########################################################################
dataSet <- "GSE34754";
chipType <- "Mapping250K_Nsp";

verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType);
print(ds);
## AffymetrixCelSet:
## Name: GSE34754
## Tags:
## Path: rawData/GSE34754/Mapping250K_Nsp
## Platform: Affymetrix
## Chip type: Mapping250K_Nsp
## Number of arrays: 12
## Names: GSM854615, GSM854616, GSM854617, ..., GSM854626 [12]
## Time period: 2011-03-29 20:30:55 -- 2011-04-05 16:54:26
## Total file size: 751.17MB
## RAM: 0.02MB



##########################################################################
# Data set:
# GSE12702/
#   Mapping250K_Nsp/
#    GSM318728.CEL, ..., GSM318767.CEL [40]
#   Mapping250K_Sty/
#    GSM318773.CEL, ..., GSM318812.CEL [40]
#
# Overall design:
#  Profiles were generated on two Affymetrix array chips: 
#  Nsp and Sty, each with ~250K SNPs represented. In all,
#  20 tumors and 20 paired normal samples were profiled,
#  40 profiles in all. Each tumor was centered on its 
#  corresponding normal pair to define copy number alterations
#  (CNA) in that tumor.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE12702
##########################################################################
dataSet <- "GSE12702";
chipType <- "Mapping250K_Nsp";

verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType);
print(ds);
## AffymetrixCelSet:
## Name: GSE12702
## Tags:
## Path: rawData/GSE12702/Mapping250K_Nsp
## Platform: Affymetrix
## Chip type: Mapping250K_Nsp
## Number of arrays: 40
## Names: GSM318728, GSM318729, GSM318730, ..., GSM318767 [40]
## Time period: 2008-04-03 13:52:53 -- 2008-04-10 20:08:49
## Total file size: 2506.11MB
## RAM: 0.04MB

chipType <- "Mapping250K_Sty";
ds <- downloadGeoRawDataSet(dataSet, chipType=chipType);
print(ds);
## AffymetrixCelSet:
## Name: GSE12702
## Tags:
## Path: rawData/GSE12702/Mapping250K_Sty
## Platform: Affymetrix
## Chip type: Mapping250K_Sty
## Number of arrays: 40
## Names: GSM318773, GSM318774, GSM318775, ..., GSM318812 [40]
## Time period: 2007-10-03 16:11:06 -- 2008-02-21 15:18:27
## Total file size: 2506.12MB
## RAM: 0.04MB

verbose && exit(verbose);
