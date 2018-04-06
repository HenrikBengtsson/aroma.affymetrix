###########################################################################/**
# @RdocClass AffymetrixFileSet
#
# @title "The AffymetrixFileSet class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixFileSet object represents a set of @see "AffymetrixFile"s
#  with \emph{identical} chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "AffymetrixFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("AffymetrixFileSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    reqFileClass <- "AffymetrixFile"
    lapply(files, FUN=function(df) {
      df <- Arguments$getInstanceOf(df, reqFileClass, .name="files")
    })
  } else if (inherits(files, "AffymetrixFileSet")) {
    return(as.AffymetrixFileSet(files))
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files))
  }


  extend(AromaMicroarrayDataSet(files=files, ...), c("AffymetrixFileSet",
                                             uses("AromaPlatformInterface")))
})




###########################################################################/**
# @RdocMethod as.AffymetrixFileSet
# @alias as.AffymetrixFileSet.list
# @alias as.AffymetrixFileSet.default
#
# @title "Coerce an object to an AffymetrixFileSet object"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Other arguments passed to @see "base::list.files".}
# }
#
# \value{
#   Returns an @see "AffymetrixFileSet" object.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.AffymetrixFileSet", "AffymetrixFileSet", function(object, ...) {
  object
})

setMethodS3("as.AffymetrixFileSet", "list", function(object, ...) {
  AffymetrixFileSet(object, ...)
})

setMethodS3("as.AffymetrixFileSet", "default", function(object, ...) {
  throw("Cannot coerce object to an AffymetrixFileSet object: ", mode(object))
})




###########################################################################/**
# @RdocMethod byPath
#
# @title "Defines an AffymetrixFileSet object by searching for Affymetrix files"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to the constructor of the static
#     (calling) class.}
#  \item{fileClass}{The name of the @see "GenericDataFile" class.}
# }
#
# \value{
#   Returns an @see "AffymetrixFileSet" object.
# }
#
# \section{Reserved filenames}{
#   Note that files with names starting with a period \code{.} are not
#   searched for.  The reason for this is that such files are reserved for
#   internal use of this package.  For instance, the package store average
#   signals across CEL files in a file named as \code{.average<something>.CEL}
#   in the same directory as the CEL files, and when such a directory is
#   scanned we do not want such files to be interpreted as data.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("byPath", "AffymetrixFileSet", function(static, ..., fileClass="AffymetrixFile") {
  # WORKAROUND: Note, NextMethod() passes all arguments ever passed in the
  # sequence of calls (i.e. the calls to the generic function, and each
  # of the "next" methods before reaching this one).  It is not possible
  # grab/exclude any of them when calling NextMethod().  Because of this
  # above '...' may contain several methods that are not accepted down
  # stream.  Indeed, here '...' will be passed by byPath() for GenericDataSet
  # to newInstance(static, ...), which will generate an error unless
  # .onUnknownArgs="ignore". /HB 2012-10-18
  NextMethod("byPath", fileClass=fileClass, .onUnknownArgs="ignore")
}, static=TRUE)


setMethodS3("getDefaultFullName", "AffymetrixFileSet", function(this, parent=1L, ...) {
  NextMethod("getDefaultFullName", parent=parent)
}, protected=TRUE)
