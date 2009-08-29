library("aroma.affymetrix");
log <- Arguments$getVerbose(-50, timestamp=TRUE);


dataSetName <- "Affymetrix_2009-CytoSampleData";
chipType <- "Cytogenetics_Array";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up CEL set and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

csR <- AffymetrixCelSet$byName(dataSetName, cdf=cdf, verbose=log);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic-crosstalk calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR);
print(acc);
csC <- process(acc, verbose=log);
print(csC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot allele pairs before and after calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (what in c("input", "output")) {
  plotAllelePairs(acc, array=1, what=what, verbose=log);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level summarization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- AvgSnpPlm(csC, mergeStrands=TRUE);
print(plm);
units <- fit(plm, verbose=log);
str(units);
ces <- getChipEffectSet(plm);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allele summation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
as <- AlleleSummation(ces);
cesS <- process(as, verbose=log);
print(cesS);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation and plotting
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cns <- CbsModel(cesS);
print(cns);

ce <- ChromosomeExplorer(cns);
setZooms(ce, 2^0:5);
print(ce);
process(ce, chromosomes=c(19, 22, 23), verbose=log);
