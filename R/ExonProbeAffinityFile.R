###########################################################################/**
# @RdocClass ExonProbeAffinityFile
#
# @title "The ExonProbeAffinityFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of probe affinities in exon array
#  probe-level models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeAffinityFile".}
#   \item{mergeGroups}{Specifies if the groups (exons) are merged or not for
#      these estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("ExonProbeAffinityFile", function(..., mergeGroups=FALSE) {
  extend(ProbeAffinityFile(...), "ExonProbeAffinityFile",
    mergeGroups=mergeGroups
  )
})


setMethodS3("getCellIndices", "ExonProbeAffinityFile", function(this, ...) {
  cells <- NextMethod("getCellIndices");

  # Merge groups?
  if (this$mergeGroups) {
    cells <- applyCdfGroups(cells, cdfMergeGroups);
  }

  cells;
})

setMethodS3("setMergeGroups", "ExonProbeAffinityFile", function(this, status, ...) {
  this$mergeGroups <- status;
})


############################################################################
# HISTORY:
# 2007-02-07
# o Created (based on SnpProbeAffinityFile.R following chat with HB on
#   2007-02-07).
############################################################################
