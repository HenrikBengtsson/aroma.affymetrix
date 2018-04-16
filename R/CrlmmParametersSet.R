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
# @author "HB"
#*/###########################################################################
setConstructorS3("CrlmmParametersSet", function(...) {
  extend(AromaUnitSignalBinarySet(...), "CrlmmParametersSet")
})


setMethodS3("byName", "CrlmmParametersSet", function(static, name, tags=NULL, ..., chipType=NULL, paths="crlmmData(|,.*)/") {
  suppressWarnings({
    path <- findByName(static, name=name, tags=tags, chipType=chipType,
                                           ..., paths=paths, mustExist=TRUE)
  })

  byPath(static, path=path, ...)
}, static=TRUE)

setMethodS3("byPath", "CrlmmParametersSet", function(static, ...) {
  suppressWarnings({
    NextMethod("byPath", pattern=".*,CRLMM[.]atb$$")
  })
})


setMethodS3("findUnitsTodo", "CrlmmParametersSet", function(this, ...) {
  # Look into the chip-effect file that comes last in a lexicographic
  # order, becuase that is updated last.
  names <- getFullNames(this)
  idx <- order(names, decreasing=TRUE)[1]
  df <- this[[idx]]
  findUnitsTodo(df, ...)
})
