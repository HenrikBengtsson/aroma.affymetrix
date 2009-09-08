library("aroma.affymetrix");

log <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup of annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CDF
cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full");

# Assert that an UGP annotation data file exists
gi <- getGenomeInformation(cdf);
print(gi);

# Assert that an UFL annotation data file exists
si <- getSnpInformation(cdf);
print(si);

# Assert than an ACS (probe-sequence) annotation files
acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName("HapMap270,6.0,CEU,testSet", cdf=cdf);
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2");
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
plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE);
print(plm);

if (length(findUnitsTodo(plm)) > 0) {
   # Fit CN probes quickly (~5-10s/array + some overhead)
  units <- fitCnProbes(plm, verbose=log);
  str(units);
  # int [1:945826] 935590 935591 935592 935593 935594 935595 ...

  # Fit remaining units, i.e. SNPs (~5-10min/array)
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
# Extract (thetaA, thetaB) for copy-neutral chromosomes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- getCdf(cesN);
gi <- getGenomeInformation(cdf);
units <- getUnitsOnChromosomes(gi, 1:22);
theta <- extractTheta(cesN, units=units, drop=TRUE);
str(theta);
print(theta[1:10,,]);
