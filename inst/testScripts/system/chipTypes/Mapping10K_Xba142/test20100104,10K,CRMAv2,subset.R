##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

cdf <- AffymetrixCdfFile$byChipType("Mapping10K_Xba142");
print(cdf);
gi <- getGenomeInformation(cdf);
print(gi);
si <- getSnpInformation(cdf);
print(si);

csR <- AffymetrixCelSet$byName("GSE8605", cdf=cdf);
print(csR);

# Run CRMAv2 on a subset of the arrays.  This is possible because
# CRMAv2 is truly a single-array method.
subset <- c(2,5,3);
csR <- extract(csR, subset);
print(csR);

acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2");
print(acc);
csC <- process(acc, verbose=log);
print(csC);

# Assert that the correct arrays have been processed
stopifnot(getNames(csC) == getNames(csR));


bpn <- BasePositionNormalization(csC, target="zero");
print(bpn);
csN <- process(bpn, verbose=log);
print(csN);

# Assert that the correct arrays have been processed
stopifnot(getNames(csN) == getNames(csR));


plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE);
print(plm);
fit(plm, verbose=log);

ces <- getChipEffectSet(plm);

# Assert that the correct arrays have been processed
stopifnot(getNames(ces) == getNames(csR));


fln <- FragmentLengthNormalization(ces, target="zero");
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);

# Assert that the correct arrays have been processed
stopifnot(getNames(cesN) == getNames(csR));
