###########################################################################/**
# @set "class=AffymetrixCelFile"
# @RdocMethod normalizeQuantile
#
# @title "Normalizes the probe intensities to a target empirical distribution"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The path where to save the normalized data files.}
#   \item{xTarget}{A @numeric @vector.  The empirical distribution
#     to which all arrays should be normalized to.}
#   \item{subsetToUpdate}{The indices of the probes to be updated.
#     If @NULL, all are updated.}
#   \item{typesToUpdate}{Types of probes to be updated.  For more details,
#     see argument \code{types} of \code{identifyCells()} for the
#     @see "AffymetrixCdfFile" class.}
#   \item{...}{Additional arguments passed to \code{normalizeQuantile()}.}
#   \item{overwrite}{If @TRUE, already normalized arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the normalized @see "AffymetrixCelFile" object.
# }
#
# @author "HB, KS"
#
# \seealso{
#   @see "aroma.light::normalizeQuantile"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("normalizeQuantile", "AffymetrixCelFile", function(this, path=file.path("normQuantile", getChipType(this)), xTarget, subsetToUpdate=NULL, typesToUpdate=NULL, ..., overwrite=FALSE, skip=!overwrite, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePathname(path)
  if (identical(getPath(this), path)) {
    throw("Cannot not normalize data file. Argument 'path' refers to the same path as the path of the data file to be normalized: ", path)
  }


  cdf <- getCdf(this)
  nbrOfCells <- nbrOfCells(cdf)

  # Argument 'subsetToUpdate':
  getFraction <- (length(subsetToUpdate) == 1) &&
                               (subsetToUpdate >= 0) && (subsetToUpdate < 1)
  if (!getFraction) {
    subsetToUpdate <- Arguments$getIndices(subsetToUpdate, max=nbrOfCells)
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- basename(getPathname(this))
  filename <- gsub("[.]cel$", ".CEL", filename);  # Only output upper case!
  pathname <- Arguments$getWritablePathname(filename, path=path,
                                         mustNotExist=(!overwrite && !skip))
  pathname <- AffymetrixFile$renameToUpperCaseExt(pathname)

  # Already normalized?
  if (skip && isFile(pathname)) {
    verbose && cat(verbose, "Normalized data file already exists: ", pathname)
    # CDF inheritance
    res <- fromFile(this, pathname)
    setCdf(res, cdf)
    return(res)
  }

  # Get all probe signals
  verbose && enter(verbose, "Reading probe intensities")
  x <- getData(this, fields="intensities", ..., verbose=less(verbose,2))
  x <- x$intensities
  verbose && exit(verbose)

  # Identify the subset of probes to be updated?
  if (getFraction || !is.null(typesToUpdate)) {
    verbose && enter(verbose, "Identifying probes to be updated")
    subsetToUpdate <- identifyCells(cdf, indices=subsetToUpdate,
                                types=typesToUpdate, verbose=less(verbose))
    verbose && exit(verbose)
  }

  # Normalize intensities
  verbose && enter(verbose, "Normalizing to empirical target distribution")
  x[subsetToUpdate] <- .normalizeQuantile(x[subsetToUpdate], xTarget=xTarget)
  # Not needed anymore
  subsetToUpdate <- NULL
  verbose && exit(verbose)

  # Write normalized data to file
  verbose && enter(verbose, "Writing normalized probe signals")

  # Write to a temporary file
  # (allow rename of existing one if forced)
  isFile <- (!skip && isFile(pathname))
  pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=less(verbose,10))

  # Create CEL file to store results, if missing
  verbose && enter(verbose, "Creating CEL file for results, if missing")
  createFrom(this, filename=pathnameT, path=NULL, verbose=less(verbose))
  verbose && exit(verbose)

  verbose && enter(verbose, "Writing normalized intensities")
  .updateCel(pathnameT, intensities=x)
  verbose && exit(verbose)
  verbose && exit(verbose)

  # Rename temporary file
  popTemporaryFile(pathnameT, verbose=less(verbose,10))

  # Return new normalized data file object
  res <- fromFile(this, pathname)

  # CDF inheritance
  setCdf(res, cdf)

  ## Create checksum file
  resZ <- getChecksumFile(res)

  res
}, private=TRUE)
