###########################################################################/**
# @RdocClass QualityAssessmentSet
#
# @title "The QualityAssessmentSet class"
#
# \description{
#  @classhierarchy
#
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to constructor of @see "AffymetrixCelSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS"
#*/###########################################################################
setConstructorS3("QualityAssessmentSet", function(...) {
  extend(AffymetrixCelSet(...), "QualityAssessmentSet")
})
