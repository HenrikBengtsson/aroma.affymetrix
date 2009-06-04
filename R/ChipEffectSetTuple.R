setConstructorS3("ChipEffectSetTuple", function(csList=NULL, ..., .setClass="ChipEffectSet") {
  extend(AffymetrixCelSetTuple(csList=csList, ..., .setClass=.setClass), "ChipEffectSetTuple")
})

setMethodS3("getFullNames", "ChipEffectSetTuple", function(this, ..., exclude="chipEffects") {
  NextMethod("getFullNames", this, ..., exclude=exclude);
})


##############################################################################
# HISTORY:
# 2007-03-19
# o Created.
##############################################################################
