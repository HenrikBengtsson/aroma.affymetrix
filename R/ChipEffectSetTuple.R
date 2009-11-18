setConstructorS3("ChipEffectSetTuple", function(csList=NULL, ..., .setClass="ChipEffectSet") {
  # Nothing do to?
  if (inherits(csList, "ChipEffectSetTuple")) {
    return(csList);
  }

  extend(AffymetrixCelSetTuple(csList=csList, ..., .setClass=.setClass), "ChipEffectSetTuple");
})


setMethodS3("getFullNames", "ChipEffectSetTuple", function(this, ..., exclude="chipEffects") {
  NextMethod("getFullNames", this, ..., exclude=exclude);
})


setMethodS3("as.ChipEffectSetTuple", "ChipEffectSetTuple", function(this, ...) {
  # Nothing do to
  this;
})


setMethodS3("as.ChipEffectSetTuple", "default", function(this, ...) {
  ChipEffectSetTuple(this, ...);
})


##############################################################################
# HISTORY:
# 2009-11-18
# o Added as.ChipEffectSetTuple().
# 2007-03-19
# o Created.
##############################################################################
