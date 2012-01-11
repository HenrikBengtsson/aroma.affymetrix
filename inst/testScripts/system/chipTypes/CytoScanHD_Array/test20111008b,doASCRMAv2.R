library("aroma.affymetrix");

verbose <- Arguments$getVerbose(-50, timestamp=TRUE);

dataSetName <- "Affymetrix_2011-CytoScanHD";
chipType <- "CytoScanHD_Array";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

csR <- AffymetrixCelSet$byName(dataSetName, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allele-specific CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- doASCRMAv2(csR, verbose=verbose);
print(dsList);
