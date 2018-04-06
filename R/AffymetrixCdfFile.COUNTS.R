setMethodS3("nbrOfCells", "AffymetrixCdfFile", function(this, ...) {
  as.integer(prod(getDimension(this, ...)))
})

setMethodS3("nbrOfUnits", "AffymetrixCdfFile", function(this, ...) {
  getHeader(this)$probesets
})

setMethodS3("nbrOfQcUnits", "AffymetrixCdfFile", function(this, ...) {
  getHeader(this)$qcprobesets
})


###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod nbrOfGroupsPerUnit
#
# @title "Gets the number of groups in each unit"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units of interest. If @NULL, all units are considered.}
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @vector of @integers.
# }
#
# \details{
#   Once read from file, this information is cached in memory for efficiency.
#   The cache can be cleared by calling \code{gc(cdf)}.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("nbrOfGroupsPerUnit", "AffymetrixCdfFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  sizes <- this$.nbrOfGroupsPerUnit

  if (force || is.null(sizes)) {
    # Check in file cache
    chipType <- getChipType(this)
    key <- list(method="nbrOfGroupsPerUnit", class=class(this)[1],
                chipType=chipType)
    if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
      key <- getCacheKey(this, method="nbrOfGroupsPerUnit", chipType=chipType)
    }
    dirs <- c("aroma.affymetrix", chipType)
    if (force) {
      sizes <- NULL
    } else {
      sizes <- loadCache(key=key, dirs=dirs)
    }

    if (is.null(sizes)) {
      verbose && enter(verbose, "Reading number of groups for *all* units")
      sizes <- .readCdfGroupNames(getPathname(this))
      sizes <- restruct(this, sizes, verbose=less(verbose, 5))
      sizes <- lapply(sizes, FUN=length)
      sizes <- unlist(sizes, use.names=FALSE)
      saveCache(sizes, key=key, dirs=dirs)
      verbose && exit(verbose)
    }

    this$.nbrOfGroupsPerUnit <- sizes
  }

  if (!is.null(units))
    sizes <- sizes[units]

  sizes
}, private=TRUE)


setMethodS3("nbrOfCellsPerUnitGroup", "AffymetrixCdfFile", function(this, units=NULL, ..., useNames=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  chipType <- getChipType(this)
  verbose && enter(verbose, "Getting the number of cells per unit group")
  verbose && cat(verbose, "Chip type: ", chipType)

  key <- list(method="nbrOfCellsPerUnitGroup", class=class(this)[1],
              chipType=chipType, useNames=useNames)
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="nbrOfCellsPerUnitGroup", chipType=chipType, useNames=useNames)
  }
  dirs <- c("aroma.affymetrix", chipType)
  if (force) {
    counts <- NULL
  } else {
    verbose && enter(verbose, "Checking for cached results")
    counts <- loadCache(key=key, dirs=dirs)
    if (!is.null(counts))
      verbose && cat(verbose, "Cached results found.")
    verbose && exit(verbose)
  }

  if (is.null(counts)) {
    verbose && enter(verbose, "Reading cell counts from CDF")
    counts <- .readCdfNbrOfCellsPerUnitGroup(getPathname(this))
    verbose && exit(verbose)

    counts <- restruct(this, counts, verbose=less(verbose, 5))

    if (!useNames) {
      names(counts) <- NULL
      counts <- lapply(counts, FUN=function(x) { names(x) <- NULL; x })
    }

    saveCache(counts, key=key, dirs=dirs)
  }

  # Subset?
  if (!is.null(units)) {
    counts <- counts[units]
  }

  verbose && exit(verbose)

  counts
}, private=TRUE)



setMethodS3("nbrOfCellsPerUnit", "AffymetrixCdfFile", function(this, units=NULL, ..., useNames=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  chipType <- getChipType(this)
  verbose && enter(verbose, "Getting the number of cells per unit")
  verbose && cat(verbose, "Chip type: ", chipType)

  key <- list(method="nbrOfCellsPerUnit", class=class(this)[1],
              chipType=chipType, useNames=useNames)
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="nbrOfCellsPerUnit", chipType=chipType, useNames=useNames)
  }
  dirs <- c("aroma.affymetrix", chipType)
  if (force) {
    counts <- NULL
  } else {
    verbose && enter(verbose, "Checking for cached results")
    counts <- loadCache(key=key, dirs=dirs)
    if (!is.null(counts))
      verbose && cat(verbose, "Cached results found.")
    verbose && exit(verbose)
  }

  if (is.null(counts)) {
    verbose && enter(verbose, "Getting number of cells per unit group")
    counts <- nbrOfCellsPerUnitGroup(this, units=NULL, useNames=useNames,
                                           ..., verbose=less(verbose, 1))

    verbose && enter(verbose, "Summing per unit")
    # Sum per unit
    counts <- sapply(counts, FUN=sum)
    verbose && exit(verbose)

    verbose && exit(verbose)

    saveCache(counts, key=key, dirs=dirs)
  }

  # Subset?
  if (!is.null(units)) {
    counts <- counts[units]
  }

  verbose && exit(verbose)

  counts
}, private=TRUE)
