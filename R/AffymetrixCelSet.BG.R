###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod calculateParametersGsb
#
# @title "Computes parameters for adjustment of specific binding"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{nbrOfPms}{The number of random PMs to use in estimation.}
#   \item{affinities}{A @numeric @vector of probe affinities.}
#   \item{seed}{An (optional) @integer specifying a temporary random seed
#     to be used during processing.  The random seed is set to its original
#     state when done.  If @NULL, it is not set.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \details{
#   This method is not constant in memory! /HB 2007-03-26
# }
#
# @author "KS"
#*/###########################################################################
setMethodS3("calculateParametersGsb", "AffymetrixCelSet", function(this, nbrOfPms=25000, affinities=NULL, seed=NULL, ..., verbose=FALSE) {
  if (is.null(affinities)) {
    throw("DEPRECATED: Argument 'affinities' to calculateParametersGsb() must not be NULL.")
  }

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed)
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  cdf <- getCdf(this)

  verbose && enter(verbose, "Computing parameters for adjustment of specific binding")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set the random seed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(seed)) {
    randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
    on.exit(randomSeed("reset"), add=TRUE)
    verbose && printf(verbose, "Random seed temporarily set (seed=%d, kind=\"L'Ecuyer-CMRG\")\n", seed)
  }

  verbose && enter(verbose, "Extracting PM indices")
  cells <- getCellIndices(cdf, useNames=FALSE, unlist=TRUE, verbose=less(verbose,2))
  pmCells <- cells[isPm(cdf, cache=FALSE)]
  # Not needed anymore
  cells <- NULL
  pmCells <- sort(pmCells)
  verbose && exit(verbose)

  narray <- length(this)

  #  set.seed(1)
  #  was present in original gcrma code; left in here to allow for consistency
  #  check between old and new versions

  # get a sorted random subset of PM to use in parameter estimation
  pmCells.random <- sample(pmCells, size=nbrOfPms)
  pmCells.random <- sort(pmCells.random)
  # Not needed anymore
  pmCells <- NULL;   # Not needed anymore

  # Garbage collect
  gc <- gc()
  verbose && print(verbose, gc)

  verbose && enter(verbose, "Extracting ", nbrOfPms, " random PM intensities across CEL set")
  # make sure we don't just sample from a single array; avoids problems
  # if we happened to choose a low quality or otherwise aberrant array
  iarray <- sample(seq_len(narray), size=nbrOfPms, replace=TRUE)

  # For each array, read the signals randomized for that array
  # Confirmed to give identical results. /HB 2007-03-26
  pathnames <- getPathnames(this)
  pm.random2 <- vector("double", length=nbrOfPms)
  for (ii in seq_len(narray)) {
    verbose && enter(verbose, sprintf("Array #%d of %d", ii, narray))
    # Cells to be read for this array
    idxs <- which(iarray == ii)
    cells <- pmCells.random[idxs]
    pm.random2[idxs] <- .readCel(pathnames[ii], indices=cells,
                       readIntensities=TRUE, readStdvs=FALSE)$intensities
    # Not needed anymore
    idxs <- cells <- NULL
    verbose && exit(verbose)
  } # for (ii ...)
  # Not needed anymore
  iarray <- pathnames <- NULL

  # Garbage collect
  gc <- gc()
  verbose && print(verbose, gc)

  verbose && exit(verbose)

  verbose && enter(verbose, "Extracting probe affinities and fitting linear model")

  aff <- affinities[pmCells.random]
  # Not needed anymore
  pmCells.random <- NULL;  # Not needed anymore

  # Work on the log2 scale
  pm.random2 <- log2(pm.random2);  # Minimize memory usage.
  # Garbage collect
  gc <- gc()
  verbose && print(verbose, gc)

  verbose && enter(verbose, "Fitting the GCRMA background linear model")
  verbose && str(verbose, pm.random2)
  verbose && str(verbose, aff)
  fit1 <- lm(pm.random2 ~ aff)
  verbose && print(verbose, fit1)
  verbose && exit(verbose)

  verbose && exit(verbose)

  verbose && exit(verbose)

  fit1$coefficients
}, private=TRUE)
