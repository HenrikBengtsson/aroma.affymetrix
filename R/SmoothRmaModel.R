###########################################################################/**
# @RdocClass SmoothRmaModel
#
# @title "The SmoothRmaModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Chromosomal Smoothing Robust Multichip Analysis
#  method.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of 
#              @see "SmoothMultiarrayModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("SmoothRmaModel", function(...) {
  extend(SmoothMultiarrayModel(...), "SmoothRmaModel");
})


setMethodS3("getAsteriskTags", "SmoothRmaModel", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", this, collapse=NULL, ...);

  # Replace first tags
  tags[1] <- "SRMA";

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  tags;
}, protected=TRUE)

setMethodS3("getRootPath", "SmoothRmaModel", function(this, ...) {
  "plmData";
})


setMethodS3("getFitUnitGroupFunction", "SmoothRmaModel", function(this, ...) {
  smoothWRMA;
})


##############################################################################
# HISTORY:
# 2007-09-20
# o Created.
##############################################################################
