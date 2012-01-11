library("aroma.affymetrix");

log <- Arguments$getVerbose(-3, timestamp=TRUE);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "Affymetrix-HeartBrain";
chipType <- "HuEx-1_0-st-v2";
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="coreR3,A20071112,EP");
print(cdf);

# Setup CEL set using the core CDF.
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Background correction and normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bc <- RmaBackgroundCorrection(csR);
print(bc);
csBC <- process(bc, verbose=log);
print(csBC);

qn <- QuantileNormalization(csBC, typesToUpdate="pm");
print(qn);
csN <- process(qn, verbose=log);
print(csN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level summarization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level summarization
plmList <- list(
  merge   = ExonRmaPlm(csN, mergeGroups=TRUE), # all exons together
  noMerge = ExonRmaPlm(csN, mergeGroups=FALSE) # each exon separately
);
print(plmList);

# Fit the PLMs
dummy <- lapply(plmList, FUN=fit, verbose=log);
rm(dummy);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Residuals and weights (will otherwise be calculated by FIRMA)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
resList <- lapply(plmList, FUN=calculateResiduals, verbose=log);
print(resList);

weightList <- lapply(plmList, FUN=calculateWeights, verbose=log);
print(weightList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FIRMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FIRMA can only be applied for PLMs with mergeGroups=TRUE
plm <- plmList$merge;
firma <- FirmaModel(plm);
print(firma);
fit(firma, verbose=log);
fScores <- getFirmaScores(firma);
print(fScores);
