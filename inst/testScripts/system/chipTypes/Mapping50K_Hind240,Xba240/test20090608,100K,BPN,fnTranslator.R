library("aroma.affymetrix");

log <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSetName <- "HapMap270,100K,CEU,testSet";
chipType <- "Mapping50K_Hind240";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert existence of probe-sequence annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acs <- AromaCellSequenceFile$byChipType(chipType);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
print(csR);

print(getFullNames(csR));
setFullNamesTranslator(csR, function(names, ...) {
  paste(names, "foo", sep=",");
})
print(getFullNames(csR));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Base-position normalization (on PM-only probes)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bpn <- BasePositionNormalization(csR, typesToFit="pm", tags="*,pm");
print(bpn);

csN <- process(bpn, verbose=log);
print(csN);
