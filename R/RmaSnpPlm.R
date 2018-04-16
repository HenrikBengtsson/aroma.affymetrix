###########################################################################/**
# @RdocClass RmaSnpPlm
#
# @title "The RmaSnpPlm class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RmaPlm".}
#   \item{mergeStrands}{If @TRUE, the sense and the anti-sense strands are
#      fitted together, otherwise separately.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("RmaSnpPlm", function(..., mergeStrands=FALSE) {
  extend(RmaPlm(...), c("RmaSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})


setMethodS3("getAsteriskTags", "RmaSnpPlm", function(this, collapse=NULL, ...) {
  # Returns 'RMA[,<flavor>]'
  tags <- NextMethod("getAsteriskTags", collapse=NULL)

  # Add class specific parameter tags
  if (!this$mergeStrands)
    tags <- c(tags, "+-")

  # Collapse
  tags <- paste(tags, collapse=collapse)

  tags
}, protected=TRUE)
