##########################################################################
# Allele-specific CRMAv2
##########################################################################
future::plan("multisession")
library("aroma.affymetrix")
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE19539"
chipType <- "GenomeWideSNP_6,Full"

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType)
print(csR)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsNList <- doASCRMAv2(csR, verbose=verbose)
print(dsNList)

dsN <- exportAromaUnitPscnBinarySet(dsNList)
print(dsN)
