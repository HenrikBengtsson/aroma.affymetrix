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
#   \item{...}{Arguments passed to @see "aroma.core::GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AffymetrixFileSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    reqFileClass <- "AffymetrixFile";
    base::lapply(files, FUN=function(df) {
      if (!inherits(df, reqFileClass))
        throw("Argument 'files' contains a non-", reqFileClass, 
                                                   " object: ", class(df)[1]);
    })
  } else if (inherits(files, "AffymetrixFileSet")) {
    return(as.AffymetrixFileSet(files));
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }


  extend(AromaMicroarrayDataSet(files=files, ...), c("AffymetrixFileSet",
                                             uses("AromaPlatformInterface")));
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
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.AffymetrixFileSet", "AffymetrixFileSet", function(object, ...) {
  object;
})

setMethodS3("as.AffymetrixFileSet", "list", function(object, ...) {
  AffymetrixFileSet(object, ...);
})

setMethodS3("as.AffymetrixFileSet", "default", function(object, ...) {
  throw("Cannot coerce object to an AffymetrixFileSet object: ", mode(object));
})




###########################################################################/**
# @RdocMethod fromFiles
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
#  \item{path}{The path where to search for Affymetrix files.}
#  \item{pattern}{The filename pattern for match files. 
#     If @NULL, filename extensions corresponding to known subclasses
#     of the abstract @see "AffymetrixFile" class are search for.}
#  \item{recursive}{If @TRUE, subdirectories are search recursively,
#     otherwise not.}
#  \item{...}{Optional arguments passed to the constructor of the
#     static (calling) class.}
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
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromFiles", "AffymetrixFileSet", function(static, ..., fileClass="AffymetrixFile") {
  # NextMethod() does not work here
  fromFiles.GenericDataFileSet(static, ..., fileClass=fileClass);
}, static=TRUE)




############################################################################
# HISTORY:
# 2008-05-18
# o Now "provides" AffymetrixPlatform methods.
# 2008-05-09
# o Removed clearCache() since it was identical to the inherited one.
# o BUG FIX: clearCache() called getCdf() which did not exist.
# 2007-09-14
# o AffymetrixFileSet now inherits from GenericDataFileSet.
# 2007-03-06
# o Added indexOf().
# 2007-02-15
# o Added getFullNames().
# 2007-02-07
# o Added sapply().
# 2007-01-14
# o Added a test for "unknown" (=unused) arguments to constructor.
# 2006-01-07
# o Added hasTags() and hasTag().
# o Now arguments '...' in fromFiles() are passed to the constructor of
#   the static class, i.e. via newInstance(static, ...).  Requested by KS.
# o Added argument 'alias' to constructor.
# o Added getAlias() and setAlias(), where the latter replaces setName().
# 2006-12-02
# o Added reorder().
# 2006-12-01
# o Now lapply() add data file names to returned list.
# 2006-11-20
# o Added support to override name of file set.
# o Added support for optional tags.
# 2006-11-02
# o Added getFullName(), getTags() and redefined getName().
# 2006-10-22
# o Now 'recursive' of fromFiles() defaults to FALSE.
# o Added getFiles() again.
# 2006-09-11
# o Added getPathnames().
# 2006-08-27
# o Added getFile() and getFiles().
# o Made filenames starting with a period reserved for internal use.
# 2006-08-26
# o Now getName() of a file set is inferred from the pathname:
#     path/to/<name>/chip_files/<"chip type">/
# 2006-08-21
# o Renamed 'array' to 'file'.
# o Extracted from AffymetrixCelSet.R.
# 2006-08-11
# o Added clearCache() which also clears the cache of all data file object.
# 2006-05-16
# o Redefined "[" to extract arrays.
# 2006-04-13
# o Added Rdoc comments for all methods.
# 2006-04-09
# o Now the read map is loaded automatically when fromFiles() used.
# 2006-03-30
# o Updated to new aroma.apd.
# 2006-03-18
# o Added argument 'subset' to calcAvgCellSignals() & normalizeQuantile().
# 2006-03-15
# o Now nbrOfCells() returns the number of cells for the first file only.
# o Now the fromFiles(static, ...) creates an object of the same class as 
#   the static object.
# 2006-03-04
# o Added mapping functions.
# o Added writeApd().
# 2006-03-03
# o Added lapply().
# 2006-03-02
# o Updated to deal with AffymetrixDataFile object instead of CEL directly.
# 2006-02-21
# o Letting readCelUnits() transform signals improves speed substantially.
# o Making use of new multi-array readCelUnits().
# 2006-02-20
# o Created.
############################################################################
