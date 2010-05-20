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

# Work with a single array
csR <- extract(csR, 1);
print(csR);

acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2");
print(acc);
csC <- process(acc, verbose=log);
print(csC);

bpn <- BasePositionNormalization(csC, target="zero");
print(bpn);
csN <- process(bpn, verbose=log);
print(csN);

plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE);
print(plm);
fit(plm, verbose=log);

ces <- getChipEffectSet(plm);
fln <- FragmentLengthNormalization(ces, target="zero");
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);
