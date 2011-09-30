if (interactive()) savehistory();
library("aroma.affymetrix");

verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

dataSet <- "Affymetrix-Tissues";
chipType <- "MoGene-1_0-st-v1";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="r3");
print(cdf);

csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
csR <- extract(csR, 1:6);
print(csR);
