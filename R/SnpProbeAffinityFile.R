###########################################################################/**
# @RdocClass SnpProbeAffinityFile
#
# @title "The SnpProbeAffinityFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of probe affinities in SNP probe-level models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeAffinityFile".}
#   \item{mergeStrands}{Specifies if the strands are merged or not for these
#      estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("SnpProbeAffinityFile", function(..., mergeStrands=FALSE) {
  this <- extend(ProbeAffinityFile(...), "SnpProbeAffinityFile",
    mergeStrands=mergeStrands
  );

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})


setMethodS3("getCellIndices", "SnpProbeAffinityFile", function(this, ...) {
  cells <- NextMethod("getCellIndices");

  # Merge strands?
  if (this$mergeStrands) {
    cells <- applyCdfGroups(cells, cdfMergeStrands);
  }

  cells;
})

setMethodS3("setMergeStrands", "SnpProbeAffinityFile", function(this, status, ...) {
  this$mergeStrands <- status;
})


############################################################################
# HISTORY:
# 2006-09-11
# o Created.
############################################################################
