# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# [() and [[() should be used to extract files (and nothing else)
#
# 2012-11-20
# o CLEANUP: Deprecated "[" and "[[" for AffymetrixCelFile,
#   AffymetrixCelSet CnagCfhFile, and CnagCfhSet, because they
#   should be used to subset files (and not units).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("[", "AffymetrixCelFile", function(this, units=NULL, drop=FALSE) {
  .Defunct("readUnits");
  data <- readUnits(this, units=units);
  if (drop && length(data) == 1)
    data <- data[[1]];
  data;
}, protected=TRUE, deprecated=TRUE)

setMethodS3("[[", "AffymetrixCelFile", function(this, unit=NULL) {
  .Defunct("readUnits");
  this[units=unit, drop=TRUE];
}, protected=TRUE, deprecated=TRUE)


setMethodS3("[", "AffymetrixCelSet", function(this, units=NULL, ..., drop=FALSE) {
  .Defunct("readUnits");
  res <- readUnits(this, units=units, ...);
  if (drop && length(res) == 1)
    res <- res[[1]];
  res;
}, protected=TRUE, deprecated=TRUE)

setMethodS3("[[", "AffymetrixCelSet", function(this, units=NULL, ...) {
  .Defunct("readUnits");
  this[units=units, ..., drop=TRUE];
}, protected=TRUE, deprecated=TRUE)


setMethodS3("[", "CnagCfhFile", function(this, units=NULL, alleles=NULL, drop=FALSE) {
  .Defunct("readUnits");
  data <- readUnits(this, units=units);
  if (!is.null(alleles)) {
    data <- data[, alleles, drop=drop];
  } else {
    if (drop && length(data) == 1)
      data <- data[[1]];
  }
  data;
}, protected=TRUE, deprecated=TRUE)

setMethodS3("[[", "CnagCfhFile", function(this, unit=NULL) {
  .Defunct("readUnits");
  this[units=unit, drop=TRUE];
}, protected=TRUE, deprecated=TRUE)


setMethodS3("[", "CnagCfhSet", function(this, units=NULL, ..., drop=FALSE) {
  .Defunct("readUnits");
  res <- readUnits(this, units=units, ...);
  if (drop && length(res) == 1)
    res <- res[[1]];
  res;
}, protected=TRUE, deprecated=TRUE)

setMethodS3("[[", "CnagCfhSet", function(this, units=NULL, ...) {
  .Defunct("readUnits");
  this[units=units, ..., drop=TRUE];
}, protected=TRUE, deprecated=TRUE)



setMethodS3("getParameterSet", "Model", function(this, ...) {
  .Deprecated("getParameters");
  getParameters(this, ...);
}, protected=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Misc.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("bgAdjustOptical", "AffymetrixCelSet", function(this, ...) {
  .Defunct("OpticalBackgroundCorrection");
}, private=TRUE, deprecated=TRUE)

setMethodS3("bgAdjustRma", "AffymetrixCelSet", function(this, ...) {
  .Defunct("RmaBackgroundCorrection");
}, private=TRUE, deprecated=TRUE)


setMethodS3("getUnitSizes", "AffymetrixCdfFile", function(this, ...) {
  .Deprecated("nbrOfGroupsPerUnit");
  nbrOfGroupsPerUnit(this, ...);
}, protected=TRUE, deprecated=TRUE)


# 2013-04-29 [HB]
# o Made getUnique(), createUnique(), createMonoCell(), doCRMA() defunct.
# 2011-04-15 [HB]
# o DEPRECATED: getUnique() and createUnique() are deprecated.
#   Use getUniqueCdf() and createUniqueCdf() instead.
setMethodS3("getUnique", "AffymetrixCdfFile", function(this, ...) {
  .Defunct("getUniqueCdf");
}, protected=TRUE, deprecated=TRUE)


setMethodS3("createUnique", "AffymetrixCdfFile", function(this, ...) {
  .Defunct("createUniqueCdf");
}, protected=TRUE, deprecated=TRUE)


setMethodS3("getMonoCell", "AffymetrixCdfFile", function(this, ...) {
  .Deprecated("getMonocellCdf");
}, protected=TRUE, deprecated=TRUE)


setMethodS3("createMonoCell", "AffymetrixCdfFile", function(this, ...) {
  .Defunct("createMonocellCdf");
}, protected=TRUE, deprecated=TRUE)


setMethodS3("getMatrixChipEffectFiles", "CopyNumberChromosomalModel", function(...) {
  .Deprecated("getDataFileMatrix");
  getDataFileMatrix(...);
}, protected=TRUE, deprecated=TRUE)


setMethodS3("calculateResiduals", "ProbeLevelModel", function(this, ...) {
  .Deprecated("calculateResidualSet");
  calculateResidualSet(this, ...);
}, private=TRUE)


setMethodS3("calculateResiduals", "FirmaModel", function(this, ...) {
  .Deprecated("calculateResidualSet");
  calculateResidualSet(this, ...);
}, private=TRUE)


# 2011-11-19
# o Deprecated doCRMA().  Use doCRMAv1() or doCRMAv2() instead.
setMethodS3("doCRMA", "default", function(dataSet, chipTypes=NULL, ..., logName=NULL, ram=NULL, verbose=-8) {
  .Defunct("doCRMAv1");
}, protected=TRUE, deprecated=TRUE)



setMethodS3("getExpectedOutputFiles", "MatSmoothing", function(this, ...) {
  .Deprecated("getExpectedOutputFullnames");
  fullnames <- getExpectedOutputFullnames(this, ...);

  # "Dummy" filenames
  filenames <- sprintf("%s.CEL", fullnames);

  filenames;
}, protected=TRUE, deprecated=TRUE)


# 2013-04-29
# o Made PdInfo2Cdf() defunct.
# 2012-03-23
# o CLEANUP: Deprecated PdInfo2Cdf() in favor (identical) pdInfo2Cdf(),
#   because the former does not follow the Aroma naming conventions.
PdInfo2Cdf <- function(...) {
  .Defunct("pdInfo2Cdf");
} # PdInfo2Cdf()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TODO, but still used alot internally.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getData", "AffymetrixCelFile", function(this, ...) {
##  .Deprecated("readRawData");
  readRawData(this, ...);
}, protected=TRUE, deprecated=TRUE)



############################################################################
# HISTORY:
# 2013-04-29
# o CLEANUP: Made several deprecated methods defunct.
# 2012-11-20
# o Deprecated getParameterSet() in favor of old getParameters().
# o Defuncted bgAdjust(Optical|Rma)() for AffymetrixCelSet.
# 2012-10-17
# o Deprecated getMonoCell() and createMonoCell().
# 2012-10-14
# o Created 999.DEPRECATED.R.
# 2011-02-19
# o Replaced deprecated getListOfChipEffectSets() with getSets() for
#   ChromosomalModel and SmoothMultiarrayModel.
# 2009-09-05
# o CLEAN UP: Now static methods fromChipType() and fromName() of
#   AffymetrixCelSet and other classes are defunct.  Instead, use static
#   methods byChipType() and byName() instead.
# 2008-06-06
# o Removed deprecated getTargetFunction() for FragmentLengthNormalization.
############################################################################
