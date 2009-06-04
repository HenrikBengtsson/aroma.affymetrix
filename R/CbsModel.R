###########################################################################/**
# @RdocClass CbsModel
#
# @title "The CbsModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Circular Binary Segmentation (CBS) model [1]. 
# }
# 
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "ChipEffectSetTuple".}
#   \item{...}{Arguments passed to the constructor of 
#              @see "CopyNumberSegmentationModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#  @see "CopyNumberSegmentationModel".
# }
#
# \references{
#  [1] Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. 
#      \emph{Circular binary segmentation for the analysis of array-based 
#      DNA copy number data. Biostatistics 5: 557-572, 2004.}\cr
#  [2] Venkatraman, E. S. & Olshen, A. B. 
#      \emph{A faster circular binary segmentation algorithm for the 
#      analysis of array CGH data}. Bioinformatics, 2007.\cr 
# }
#*/###########################################################################
setConstructorS3("CbsModel", function(cesTuple=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(cesTuple)) {
    require("DNAcopy") || throw("Package not loaded: DNAcopy");
  }

  extend(CopyNumberSegmentationModel(cesTuple=cesTuple, ...), "CbsModel")
})

setMethodS3("getFitFunction", "CbsModel", function(this, ...) {
  segmentByCBS;
});


##############################################################################
# HISTORY:
# 2009-05-16
# o Added getFitFunction().  Removed fitOne().
# 2007-08-20
# o Created from GladModel.R.
##############################################################################
