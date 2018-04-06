###########################################################################/**
# @RdocClass AffymetrixFile
#
# @title "The abstract AffymetrixFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixFile object represents a single Affymetrix file,
#  e.g. an Affymetrix CEL file or an Affymetrix CDF file.
#  Note that this class is abstract and can not be instanciated, but
#  instead you have to use one of the subclasses or the generic
#  \code{fromFile()} method.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   An object of this class is typically part of an @see "AffymetrixFileSet".
# }
#*/###########################################################################
setConstructorS3("AffymetrixFile", function(...) {
  extend(AromaMicroarrayDataFile(...), c("AffymetrixFile",
                                             uses("AromaPlatformInterface")))
}, abstract=TRUE)
