###########################################################################/**
# @RdocClass SnpPlm
#
# @title "The SnpPlm interface class"
#
# \description{
#  @classhierarchy
#
#  An @see "R.oo::Interface" implementing methods special for
#  @see "ProbeLevelModel"s specific to SNP arrays.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# \section{Requirements}{
#   Classes inheriting from this @see "R.oo::Interface" must provide the
#   following fields:
#   \itemize{
#    \item{mergeStrands}{A @logical value indicating if strands should be
#       merged or not.}
#   }
# }
#
# @examples "../incl/SnpPlm.Rex"
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SnpPlm", function(...) {
  extend(Interface(), "SnpPlm")
})


setMethodS3("getParameters", "SnpPlm", function(this, ...) {
  params <- NextMethod("getParameters")
  params$mergeStrands <- this$mergeStrands
  params
}, protected=TRUE)


setMethodS3("getCellIndices", "SnpPlm", function(this, ..., verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  cdfMergeStrands <- affxparser::cdfMergeStrands


  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)


  verbose && enter(verbose, "Identifying cell indices for a SnpPlm")

  cells <- NextMethod("getCellIndices", verbose=verbose)

  # Merge strands?
  if (this$mergeStrands) {
    verbose && enter(verbose, "Merging strands")
    cells <- .applyCdfGroups(cells, cdfMergeStrands)
    verbose && exit(verbose)
  }

  verbose && exit(verbose)

  cells
})

setMethodS3("getChipEffectSetClass", "SnpPlm", function(this, ...) {
  SnpChipEffectSet
}, private=TRUE)


setMethodS3("getChipEffectSet", "SnpPlm", function(this, ...) {
  ces <- NextMethod("getChipEffectSet")
  setMergeStrands(ces, this$mergeStrands)
  ces
})

setMethodS3("getProbeAffinityFile", "SnpPlm", function(this, ..., .class=SnpProbeAffinityFile) {
  paf <- NextMethod("getProbeAffinityFile", .class=.class)
  setMergeStrands(paf, this$mergeStrands)
  paf
})

setMethodS3("getMergeStrands", "SnpPlm", function(this, ...) {
  this$mergeStrands
})

setMethodS3("setMergeStrands", "SnpPlm", function(this, status, ...) {
  # Argument 'status':
  status <- Arguments$getLogical(status)

  oldStatus <- getCombineAlleles(this)

  ces <- getChipEffectSet(this)
  setMergeStrands(ces, status, ...)
  paf <- getProbeAffinityFile(this)
  setMergeStrands(paf, status, ...)
  this$mergeStrands <- status

  invisible(oldStatus)
})
