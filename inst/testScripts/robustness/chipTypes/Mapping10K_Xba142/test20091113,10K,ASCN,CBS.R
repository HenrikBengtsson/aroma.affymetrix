############################################################################
# ROBUSTNESS TEST
############################################################################
library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";

cdf <- AffymetrixCdfFile$byChipType(chipType);
cs <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
keep <- 1:6;
cs <- extract(cs, keep);
print(cs);

acc <- AllelicCrosstalkCalibration(cs);
print(acc);
csC <- process(acc, verbose=log);
print(csC);

plm <- RmaSnpPlm(csC, mergeStrands=TRUE, shift=300);
print(plm);
fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);

fln <- FragmentLengthNormalization(ces);
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);

# This should throw "Exception: Unsupported chip effects. ..."
cbs <- NULL;
tryCatch({
  cbs <- CbsModel(cesN);
}, error = function(ex) {
  print(ex);
})
print(cbs);
stopifnot(is.null(cbs));


############################################################################
# HISTORY:
# 2009-11-13
# o Created.  Thanks to Pierre Neuvial for the report.
############################################################################

