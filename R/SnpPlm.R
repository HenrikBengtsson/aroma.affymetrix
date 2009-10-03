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
# @author
#*/###########################################################################
setConstructorS3("SnpPlm", function(...) {
  extend(Interface(), "SnpPlm");
})

## setMethodS3("getSubname", "SnpPlm", function(this, ...) {
##   s <- NextMethod("getSubname", this, ...);
##   if (this$mergeStrands) {
##     s <- sprintf("%sStrandless", s);
##   } else {
##     s <- sprintf("%sStrands", s);
##   }
##   s;
## })


setMethodS3("getParameterSet", "SnpPlm", function(this, ...) {
  params <- NextMethod("getParameterSet", this, ...);
  params$mergeStrands <- this$mergeStrands;
  params;
}, private=TRUE)


setMethodS3("getCellIndices", "SnpPlm", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Identifying cell indices for a SnpPlm");

  cells <- NextMethod("getCellIndices", this, ..., verbose=verbose);

  # Merge strands?
  if (this$mergeStrands) {
    verbose && enter(verbose, "Merging strands");
    cells <- applyCdfGroups(cells, cdfMergeStrands);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);
  
  cells;
})

setMethodS3("getChipEffectSetClass", "SnpPlm", function(this, ...) {
  SnpChipEffectSet;
}, private=TRUE)


setMethodS3("getChipEffectSet", "SnpPlm", function(this, ...) {
  ces <- NextMethod("getChipEffectSet", this, ...);
  setMergeStrands(ces, this$mergeStrands);
  ces;
})

setMethodS3("getProbeAffinityFile", "SnpPlm", function(this, ..., .class=SnpProbeAffinityFile) {
  paf <- NextMethod("getProbeAffinityFile", this, ..., .class=.class);
  setMergeStrands(paf, this$mergeStrands);
  paf;
})

setMethodS3("getMergeStrands", "SnpPlm", function(this, ...) {
  this$mergeStrands;
})

setMethodS3("setMergeStrands", "SnpPlm", function(this, status, ...) {
  # Argument 'status':
  status <- Arguments$getLogical(status);

  oldStatus <- getCombineAlleles(this);

  ces <- getChipEffectSet(this);
  setMergeStrands(ces, status, ...);
  paf <- getProbeAffinityFile(this);
  setMergeStrands(paf, status, ...);
  this$mergeStrands <- status;

  invisible(oldStatus);
})


############################################################################
# HISTORY:
# 2009-02-03
# o BUG FIX: setMergeStrands() of SnpPlm did not update the setting of the
#   SnpPlm itself, only the underlying parameter files.
# o Added getMergeStrands() to SnpPlm.
# 2006-09-11
# o The intention is to use SnpPlm as an interface class (that is a class
#   that must not have any fields!) but any class "implementing" this class 
#   (by adding it to its list of classes), will have these methods too.
# o Created from the RmaSnpPlm.
############################################################################
