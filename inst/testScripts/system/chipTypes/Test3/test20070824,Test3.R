library("aroma.affymetrix")
log <- Verbose(threshold=-4, timestamp=TRUE);

dataSetName <- "FusionSDK_Test3";
chipType <- "Test3";

cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
print(cs);

cdf <- getCdf(cs);
print(cdf);

ae <- ArrayExplorer(cs);
setColorMaps(ae, "log2,yellow");
process(ae, verbose=log);
