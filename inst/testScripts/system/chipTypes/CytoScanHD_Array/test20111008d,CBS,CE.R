library("aroma.affymetrix");

verbose <- Arguments$getVerbose(-50, timestamp=TRUE);

dataSetName <- "Affymetrix_2011-CytoScanHD";
tags <- "ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY";
chipType <- "CytoScanHD_Array";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- AromaUnitTotalCnBinarySet$byName(dataSetName, tags=tags, chipType=chipType);
print(ds);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
seg <- CbsModel(ds);
print(seg);

ce <- ChromosomeExplorer(seg);
print(ce);

process(ce, arrays=1, chromosomes=1:2, verbose=verbose);
