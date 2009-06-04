###########################################################################/**
# @RdocClass CnProbeAffinityFile
#
# @title "The CnProbeAffinityFile class"
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
#   \item{...}{Arguments passed to @see "SnpProbeAffinityFile".}
#   \item{combineAlleles}{If @FALSE, allele A and allele B are treated 
#      seperately, otherwise together.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("CnProbeAffinityFile", function(..., combineAlleles=FALSE) {
  this <- extend(SnpProbeAffinityFile(...), "CnProbeAffinityFile",
    combineAlleles=combineAlleles
  );

  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})


setMethodS3("getCellIndices", "CnProbeAffinityFile", function(this, ...) {
  cells <- NextMethod("getCellIndices", this, ...);

  # If combining alleles, return only every second group.
  # In order to improve readability we merge the names of alleles groups
  # combined, e.g. groups 'C' and 'G' become group 'CG'.
  if (this$combineAlleles) {
    cells <- applyCdfGroups(cells, function(groups) {
      ngroups <- length(groups);
      odds <- seq(from=1, to=ngroups, by=2);
      names <- names(groups);
      groups <- groups[odds];
      if (ngroups >= 2) {
        evens <- seq(from=2, to=ngroups, by=2);
        names <- paste(names[odds], names[evens], sep="");
        names(groups) <- names;
      }
      groups;
    })
  }

  cells;
})

setMethodS3("setCombineAlleles", "CnProbeAffinityFile", function(this, status, ...) {
  this$combineAlleles <- status;
})



############################################################################
# HISTORY:
# 2006-09-12
# o Updated.
# 2006-09-11
# o Created.
############################################################################
