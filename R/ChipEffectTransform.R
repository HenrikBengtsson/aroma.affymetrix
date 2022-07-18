###########################################################################/**
# @RdocClass ChipEffectTransform
#
# @title "The ChipEffectTransform class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a transform that transforms chip-effect
#  estimates obtained from probe-level modeling.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{The input data set as an @see "ChipEffectSet".}
#   \item{...}{Arguments passed to the constructor of @see "Transform".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   Subclasses must implement the \code{process()} method.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("ChipEffectTransform", function(dataSet=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "ChipEffectSet")
  }

  extend(Transform(dataSet=dataSet, ...), "ChipEffectTransform")
}, abstract=TRUE)


setMethodS3("getRootPath", "ChipEffectTransform", function(this, ...) {
  "plmData"
}, protected=TRUE)
