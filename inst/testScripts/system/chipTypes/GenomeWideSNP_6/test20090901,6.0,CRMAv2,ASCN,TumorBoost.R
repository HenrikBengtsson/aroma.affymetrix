library("aroma.affymetrix");

log <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup of annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CDF
cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6");

# Assert that an UGP annotation data file exists
gi <- getGenomeInformation(cdf);
print(gi);

# Assert that an UFL annotation data file exists
si <- getSnpInformation(cdf);
print(si);

# Assert than an ACS (probe-sequence) annotation files
acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE));
print(acs);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "TCGA,OV,testSet,pairs,Broad";
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRMA v2
tags <- c("*");
alpha <- c(0.1, 0.075, 0.05, 0.03, 0.01, 0.0025, 0.001, 1e-04);
if (FALSE) {
  # CRMA v1
  alpha <- c(0.1, 0.075, 0.05, 0.03, 0.01);
  tags <- c(tags, sprintf("alpha=%.2f", alpha[length(alpha)]));
}
acc <- AllelicCrosstalkCalibration(csR, alpha=alpha, pairBy="sequence", tags=tags);
print(acc);
csC <- process(acc, verbose=log);
print(csC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Base-position normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bpn <- BasePositionNormalization(csC, target="zero");
print(bpn);

csN <- process(bpn, verbose=log);
print(csN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allele-specific chip effect estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- AvgSnpPlm(csN, mergeStrands=TRUE);
print(plm);

if (length(findUnitsTodo(plm)) > 0) {
  # Fit CN probes quickly (~5-10s/array + some overhead)
  units <- fitCnProbes(plm, verbose=log);
  str(units);
  units <- fit(plm, verbose=log);
  str(units);
}
ces <- getChipEffectSet(plm);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces, target="zero");
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Export (theta,beta) for all arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- exportTotalAndFracB(cesN, verbose=log);
print(dsList);
