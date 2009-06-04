library("aroma.affymetrix");

log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSet <- "Affymetrix-HeartBrain";
chipType <- "HG-U133_Plus_2";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);

bg <- GcRmaBackgroundCorrection(csR);
print(bg);
csB <- process(bg, verbose=log);
print(csB);

