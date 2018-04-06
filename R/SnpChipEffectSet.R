###########################################################################/**
# @RdocClass SnpChipEffectSet
#
# @title "The SnpChipEffectSet class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in the probe-level models.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ChipEffectSet".}
#   \item{mergeStrands}{Specifies if the strands are merged or not for these
#      estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SnpChipEffectSet", function(..., mergeStrands="byFirstFile") {
  this <- extend(ChipEffectSet(...), "SnpChipEffectSet")
  setMergeStrands(this, mergeStrands)
  this
})


setMethodS3("byPath", "SnpChipEffectSet", function(static, ..., mergeStrands="auto") {
  NextMethod("byPath", mergeStrands=mergeStrands)
}, static=TRUE, protected=TRUE)



setMethodS3("getAverageFile", "SnpChipEffectSet", function(this, ...) {
  res <- NextMethod("getAverageFile")
  res$mergeStrands <- getMergeStrands(this)
  res
})



setMethodS3("getChipEffectFileClass", "SnpChipEffectSet", function(static, ...) {
  SnpChipEffectFile
}, static=TRUE, private=TRUE)

setMethodS3("getMergeStrands", "SnpChipEffectSet", function(this, ...) {
  if (length(this) == 0L)
    return(FALSE)
  ce <- getOneFile(this)
  ce$mergeStrands
})

setMethodS3("setMergeStrands", "SnpChipEffectSet", function(this, status, ...) {
  if (length(this) == 0)
    return(FALSE)

  oldStatus <- getMergeStrands(this)

  if (identical(status, "byFirstFile")) {
    status <- oldStatus
  }

  if (identical(status, "auto")) {
    status <- inferParameters(this, ...)$mergeStrands
  }

  # Argument 'status':
  status <- Arguments$getLogical(status)

  # Update all chip-effect files
  lapply(this, FUN=function(ce) {
    ce$mergeStrands <- status
  })

  invisible(oldStatus)
})


setMethodS3("inferParameters", "SnpChipEffectSet", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Infer (mergeStrands) parameters from stored data in quartet units in CEL set")

  # Identify units with quartets
  cdf <- getCdf(this)
  cdfPathname <- getPathname(cdf)
  nbrOfUnits <- .readCdfHeader(cdfPathname)$nunits
  allUnits <- seq_len(nbrOfUnits)

  ce <- getOneFile(this, mustExist=TRUE)
  cePathname <- getPathname(ce)

  verbose && cat(verbose, "Pathname: ", cePathname)

  mergeStrands <- NA
  while(length(allUnits) > 0L) {
    verbose && cat(verbose, "Units left: ", length(allUnits))
    uu <- seq_len(min(10e3,length(allUnits)))
    units <- allUnits[uu]
    allUnits <- allUnits[-uu]
    # Not needed anymore
    uu <- NULL

    # Identify units that are quartets
    unitSizes <- .readCdfGroupNames(cdfPathname, units=units)
    names(unitSizes) <- NULL
    unitSizes <- sapply(unitSizes, FUN=length)
    units <- units[unitSizes == 4L]
    # Not needed anymore
    unitSizes <- NULL
    if (length(units) > 0L) {
      verbose && cat(verbose, "Scanning units:")
      verbose && str(verbose, units)
      # Infer parameters from 'intensities'
      values <- .readCelUnits(cePathname, units=units,
                   readIntensities=TRUE, readStdvs=FALSE, dropArrayDim=TRUE)
      # Not needed anymore
      units <- NULL

      # Put quartets by columns
      values <- matrix(unlist(values, use.names=FALSE), nrow=4L)

      # Keep only estimated units without NAs
      csums <- colSums(values)
      values <- values[,is.finite(csums) & (csums > 0),drop=FALSE]
      # Not needed anymore
      csums <- NULL
      verbose && cat(verbose, "Values quartets:")
      verbose && print(verbose, values[,seq_len(min(ncol(values),6)),drop=FALSE])
      if (ncol(values) > 0L) {
        t <- .rowMeans(values)
        if (length(t) > 0L) {
          isZero <- isZero(t)
          if (!all(isZero)) {
            mergeStrands <- isZero[3] && isZero[4]
            break
          }
          # Not needed anymore
          isZero <- NULL
        }
        # Not needed anymore
        t <- NULL
      }
      # Not needed anymore
      values <- NULL
    }
  } # while(...)

  if (is.na(mergeStrands)) {
    throw("Failed to infer parameter 'mergeStrands' from chip-effect file: ", cePathname)
  }

  res <- list(mergeStrands=mergeStrands)

  verbose && str(verbose, res)
  verbose && exit(verbose)

  res
}, private=TRUE)
