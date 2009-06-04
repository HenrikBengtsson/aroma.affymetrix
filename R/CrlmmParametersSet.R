###########################################################################/**
# @RdocClass CrlmmParametersSet
#
# @title "The CrlmmParametersSet class"
#
# \description{
#  @classhierarchy
#
#  An CrlmmParametersSet object represents a set of 
#  @see "CrlmmParametersFile"s with \emph{identical} chip types.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to 
#     @see "aroma.core::AromaUnitSignalBinarySet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("CrlmmParametersSet", function(...) {
  extend(AromaUnitSignalBinarySet(...), "CrlmmParametersSet");
})


setMethodS3("byName", "CrlmmParametersSet", function(static, name, tags=NULL, ..., chipType=NULL, paths="crlmmData") {
  suppressWarnings({
    path <- findByName(static, name=name, tags=tags, chipType=chipType, 
                                           ..., paths=paths, mustExist=TRUE);
  })

  fromFiles(static, path=path, ...);
}, static=TRUE) 

setMethodS3("fromFiles", "CrlmmParametersSet", function(static, ...) {
  suppressWarnings({
    fromFiles.GenericDataFileSet(static, ..., pattern=".*,CRLMM[.]atb$$");
  })
})


setMethodS3("findUnitsTodo", "CrlmmParametersSet", function(this, ...) {
  # Look into the last chip-effect file since that is updated last
  ce <- getFile(this, length(this));
  findUnitsTodo(ce, ...);
})



############################################################################
# HISTORY:
# 2008-12-08
# o Added findUnitsTodo() and extractCalls().
# 2008-12-06
# o Created.
############################################################################
