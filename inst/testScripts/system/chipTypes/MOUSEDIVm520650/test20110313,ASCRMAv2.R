library("aroma.affymetrix")

verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

dataSetName <- "Affymetrix-MouseDiversity,Set1of5,testSet";
chipType <- "MOUSEDIVm520650";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up CEL set and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

ugp <- getAromaUgpFile(cdf);
print(ugp);

ufl <- getAromaUflFile(cdf);
print(ufl);

acs <- getAromaCellSequenceFile(cdf);
print(acs);

csR <- AffymetrixCelSet$byName(dataSetName, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- doASCRMAv2(csR, lengthRange=c(450,2000), verbose=verbose);
print(dsList);
