###########################################################################/**
# @RdocClass ExonChipEffectSet
#
# @title "The ExonChipEffectSet class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in the probe-level models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ChipEffectSet".}
#   \item{mergeGroups}{Specifies if groups (individual exons in a CDF
#        file) are merged or not for these estimates, i.e. whether
#        transcript-level expression is to be estimated.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("ExonChipEffectSet", function(..., mergeGroups=TRUE) {
  this <- extend(ChipEffectSet(...), "ExonChipEffectSet");
  setMergeGroups(this, mergeGroups);
  this;
})

setMethodS3("byPath", "ExonChipEffectSet", function(static, ..., mergeGroups="auto") {
  byPath.ChipEffectSet(static, ..., mergeGroups=mergeGroups);
}, protected=TRUE, static=TRUE)



setMethodS3("getAverageFile", "ExonChipEffectSet", function(this, ...) {
  res <- NextMethod(generic="getAverageFile", object=this, ...);
  res$mergeGroups <- getMergeGroups(this);
  res;
})



setMethodS3("getChipEffectFileClass", "ExonChipEffectSet", function(static, ...) {
  ExonChipEffectFile;
}, static=TRUE, private=TRUE)

setMethodS3("getMergeGroups", "ExonChipEffectSet", function(this, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);
  ce <- getFile(this, 1);
  ce$mergeGroups;
})

setMethodS3("setMergeGroups", "ExonChipEffectSet", function(this, status, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);

  oldStatus <- getMergeGroups(this);

#  if (identical(status, "auto"))
#    status <- inferParameters(this)$mergeGroups;

  # Argument 'status':
  status <- Arguments$getLogical(status);

  # Update all chip-effect files
  lapply(this, function(ce) {
    ce$mergeGroups <- status;
  })

  invisible(oldStatus);
})

setMethodS3("getFirstCellPerUnitIndices", "ExonChipEffectSet", function(this, ...) {

  cdf <- getCdf(this);
  idx <- getFirstCellIndices(cdf, ...);
  idx <- base::lapply(base::lapply(idx, .subset2, 1), .subset2, 1);
  idx <- unlist(idx, use.names=FALSE);
  return(idx);
  
})


setMethodS3("findUnitsTodo", "ExonChipEffectSet", function(this, ...) {
  # Look into the last chip-effect file since that is updated last
  ece <- getFile(this, length(this));
  findUnitsTodo(ece, ...);
})



############################################################################
# HISTORY:
# 2008-05-08
# o Made fromFiles() protected.
# 2007-02-08
# o Created (based on SnpChipEffectSet.R following chat with HB on
#   2007-02-07).
############################################################################
