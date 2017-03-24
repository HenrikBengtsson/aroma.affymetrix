# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFUNCT
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEPRECATED
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TODO, but still used alot internally.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getData", "AffymetrixCelFile", function(this, ...) {
##  .Deprecated("readRawData");
  readRawData(this, ...);
}, protected=TRUE, deprecated=TRUE)


setMethodS3("nbrOfArrays", "AffymetrixCelSet", function(this, ...) {
##  .Deprecated("length")
  length(this, ...)
}, protected=TRUE)

setMethodS3("nbrOfArrays", "AffymetrixCnChpSet", function(this, ...) {
##  .Deprecated("length")
  length(this, ...);
}, protected=TRUE)

setMethodS3("nbrOfArrays", "CnagCfhSet", function(this, ...) {
##  .Deprecated("length")
  length(this, ...);
}, protected=TRUE)

setMethodS3("nbrOfArrays", "DChipDcpSet", function(this, ...) {
##  .Deprecated("length")
  length(this, ...);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2015-12-04
# o Moved more nbrOfArrays() methods to here. Still to be deprecated.
# 2014-02-28
# o Remove previously defunct methods.
# o CLEANUP: Defuncted deprecated patch() for AromaAffymetrix and
#   getMonoCell() & getUnitSizes() for AffymetrixCdfFile.
# 2013-10-07
# o CLEANUP: Deprecated patch() for AromaAffymetrix.
# 2014-01-04
# o Dropped defunct '[' and '[[' for AffymetrixCelSet and CnagCfhSet,
#   making ditto for super classes to be used.
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
