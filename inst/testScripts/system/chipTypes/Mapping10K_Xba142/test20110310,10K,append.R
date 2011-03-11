library("aroma.affymetrix");

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);
csR <- append(csR, csR);
print(csR);


