# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: Platform specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("importFrom", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  if (inherits(src, "AffymetrixNetAffxCsvFile")) {
    importFromAffymetrixNetAffxCsvFile(this, src, ...)
  } else if (inherits(src, "DChipGenomeInformation")) {
    importFromDChipGenomeInformation(this, src, ...)
  } else if (inherits(src, "GenomeInformation")) {
    importFromGenomeInformation(this, src, ...)
  } else if (inherits(src, "AffymetrixTabularFile")) {
    importFromAffymetrixTabularFile(this, src, ...)
  } else if (inherits(src, "GenericTabularFile")) {
    importFromGenericTabularFile(this, src, ...)
  } else {
    throw("Do not know how to import from an src of class ", class(src)[1])
  }
})



setMethodS3("getCdf", "AromaUnitTabularBinaryFile", function(this, ..., force=FALSE, .old=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  cdf <- this$.cdf
  if (force || is.null(cdf)) {
    if (.old) {
      # Generate all possible fullname 'chipTypes' and search for the existance
      # of a CDF with the longest name *and* that have the same number of units.
      # This is AD HOC! We should really store the full chiptype of the
      # corresponding CDF in the internal file header. /HB 2007-12-10

      verbose && enter(verbose, "Searching for a match CDF")

      verbose && cat(verbose, "Filename: ", getFilename(this))
      name <- getName(this, ...)
      tags <- getTags(this, collapse=NULL, ...)

      validator <- function(cdf, ...) {
        (nbrOfUnits(cdf) == nbrOfUnits(this))
      }
      pathname <- findByCdf2(chipType=name, tags=tags, validator=validator,
                                                    verbose=less(verbose, 1))
      if (is.null(pathname)) {
        throw("Failed to locate a CDF for ", class(this)[1],
              " that have ", nbrOfUnits, " units: ", getFullName(this))
      }

      cdf <- AffymetrixCdfFile$fromFile(pathname)

      verbose && exit(verbose)
    } else {
      chipType <- getChipType(this)
      nbrOfUnits <- nbrOfUnits(this)
      cdf <- AffymetrixCdfFile$byChipType(chipType, nbrOfUnits=nbrOfUnits)
    }

    this$.cdf <- cdf
  }

  cdf
})



###########################################################################/**
# @set "class=AromaUnitTabularBinaryFile"
# @RdocMethod allocateFromCdf
#
# @title "Creates an AromaUnitTabularBinaryFile mapping to a given CDF"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{The @see "AffymetrixCdfFile" used as a template and from
#      which the (full) chip type is taken.}
#   \item{...}{Additional arguments passed to \code{allocate()} of
#      @see "aroma.core::AromaTabularBinaryFile".}
# }
#
# \value{
#  Returns a @see "aroma.core::AromaUnitTabularBinaryFile" object.
# }
#
# @author "HB"
#
# \seealso{
#   To update to file footer afterwards, see \code{writeFooter()}.
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("allocateFromCdf", "AromaUnitTabularBinaryFile", function(static, cdf, ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile")

  allocateFromUnitNamesFile(static, unf=cdf, ...)
}, static=TRUE)




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: Platform specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


setMethodS3("importFromGenericTabularFile", "AromaUnitTabularBinaryFile", abstract=TRUE)


setMethodS3("importFromAffymetrixTabularFile", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  # Argument 'src':
  src <- Arguments$getInstanceOf(src, "AffymetrixTabularFile")

  importFromGenomeInformation(this, src, ...)
})


setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE)

setMethodS3("importFromDChipGenomeInformation", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  # Argument 'src':
  src <- Arguments$getInstanceOf(src, "DChipGenomeInformation")

  importFromGenomeInformation(this, src, ...)
})


setMethodS3("importFromGenomeInformation", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE)
