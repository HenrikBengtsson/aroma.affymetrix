###########################################################################/**
# @RdocClass ChipEffectSet
#
# @title "The ChipEffectSet class"
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
#   \item{...}{Arguments passed to @see "AffymetrixCelSet".}
#   \item{probeModel}{The specific type of model, e.g. \code{"pm"}.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   An object of this class is typically obtained through the
#   \code{getChipEffectSet()} method for the @see "ProbeLevelModel" class.
# }
#
#*/###########################################################################
setConstructorS3("ChipEffectSet", function(..., probeModel=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probeModel':
  probeModel <- match.arg(probeModel)

  extend(ParameterCelSet(...), "ChipEffectSet",
    "cached:.firstCells" = NULL,
    probeModel = probeModel
  )
})


setMethodS3("as.character", "ChipEffectSet", function(x, ...) {
  # To please R CMD check
  this <- x

  s <- NextMethod("as.character")
  s <- c(s, sprintf("Parameters: %s", getParametersAsString(this)))
  s
}, protected=TRUE)



setMethodS3("getChipEffectFileClass", "ChipEffectSet", function(static, ...) {
  ChipEffectFile
}, static=TRUE, private=TRUE)


setMethodS3("findByName", "ChipEffectSet", function(static, ..., paths="plmData(|,.*)/") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'paths':
  if (is.null(paths)) {
    paths <- eval(formals(findByName.ChipEffectSet)[["paths"]], enclos = baseenv())
  }

  NextMethod("findByName", paths=paths)
}, static=TRUE, protected=TRUE)


setMethodS3("byPath", "ChipEffectSet", function(static, path="plmData/", pattern=",chipEffects[.](c|C)(e|E)(l|L)(|[.]lnk|[.]LNK)$", cdf=NULL, checkChipType=FALSE, ..., fileClass=NULL) {
  # Argument 'cdf':
  if (!is.null(cdf)) {
    cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile")
    if (!isMonocellCdf(cdf)) {
      cdf <- getMonocellCdf(cdf)
    }
  }

  # Argument 'fileClass':
  if (is.null(fileClass)) {
    fileClass <- gsub("Set$", "File", class(static)[1])
  }

  ## Don't explicitly pass the first argument after 'static', otherwise
  ## it (here argument 'path') will be part of '...' as well. /HB 2013-07-28
  ## Example: ChipEffectSet$byPath(path, private=TRUE) gave "## Error in
  ## list.files(path = path, pattern = pattern, full.names = TRUE,  :
  ## invalid 'recursive' argument", because 'recursive' was also assigned
  ## the value of 'path' due to a bug(?) in how NextMethod() handles
  ## '...', cf. R-devel thread "Do *not* pass '...' to NextMethod() - it'll
  ## do it for you; missing documentation, a bug or just me?" on 2012-10-16
  ## [https://stat.ethz.ch/pipermail/r-devel/2012-October/065016.html]
  ## NB: This is actually the same fix as done in AffymetrixCnChpSet$byPath().
  NextMethod("byPath", pattern=pattern, fileClass=fileClass, cdf=cdf, checkChipType=checkChipType)
}, static=TRUE, protected=TRUE)



setMethodS3("fromDataSet", "ChipEffectSet", function(static, dataSet, path, name=getName(dataSet), cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  # Get the ChipEffectFile class specific for this set
  clazz <- getChipEffectFileClass(static)

  verbose && enter(verbose, "Retrieving chip-effects from data set")
  ces <- vector("list", length(dataSet))
  verbose && cat(verbose, "Data set: ", name)
  for (kk in seq_along(dataSet)) {
    df <- dataSet[[kk]]
    verbose && enter(verbose, sprintf("Retrieving chip-effect #%d of %d (%s)",
                                               kk, length(ces), getName(df)))
    ce <- clazz$fromDataFile(df, path=path, name=name, cdf=cdf, ...,
                                                       verbose=less(verbose))
    if (is.null(cdf)) {
      verbose && enter(verbose, "Retrieving the CDF for the chip-effect file")
      cdf <- getCdf(ce)
      verbose && exit(verbose)
    }
    ces[[kk]] <- ce
    verbose && exit(verbose)
  }
  verbose && exit(verbose)

  # Create an ChipEffectSet
  newInstance(static, ces)
}, static=TRUE, protected=TRUE)


setMethodS3("getCellIndices", "ChipEffectSet", function(this, ...) {
  # Use the first chip-effect file to get the CDF structure.
  # Note: Ideally we want to define a special CDF class doing this
  # instead of letting the data file do this. /HB 2006-12-18
  ce <- getOneFile(this)
  getCellIndices(ce, ...)
})


setMethodS3("readUnits", "ChipEffectSet", function(this, units=NULL, cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Reading chip effects unit by unit for ", length(this), " arrays")

  if (is.null(cdf)) {
    verbose && enter(verbose, "Getting cell indices from CDF")
    cdf <- getCellIndices(this, units=units, verbose=less(verbose))
    verbose && exit(verbose)
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  verbose && enter(verbose, "Calling readUnits() in superclass")
  res <- NextMethod("readUnits", units=cdf, verbose=less(verbose))
  verbose && exit(verbose)

  # Get first chip-effect file and use that to decode the read structure
  # This takes some time for a large number of units /HB 2006-10-04
  ce <- getOneFile(this)
  res <- decode(ce, res, verbose=less(verbose))

  verbose && exit(verbose)

  res
})


setMethodS3("updateUnits", "ChipEffectSet", function(this, units=NULL, cdf=NULL, data, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Updating chip-effect files")

  # Get the CDF structure for all chip-effect files
  if (is.null(cdf)) {
    cdf <- getCellIndices(this, units=units, verbose=less(verbose, 1))
  }

  # Update each file one by one
  arrays <- seq_along(this)
  nbrOfArrays <- length(arrays)
  verbose && cat(verbose, "Number of files: ", nbrOfArrays)

  if (nbrOfArrays > 1L) {
    verbose && enter(verbose, "Making sure the files are updated in lexicographic order")
    # Reorder such that the file with the "last" name is saved last
    fullnames <- getFullNames(this)
    o <- order(fullnames, decreasing=FALSE)
    arrays <- arrays[o]
    verbose && str(verbose, arrays)
    verbose && cat(verbose, "Last array: ", fullnames[arrays[nbrOfArrays]])
    # Not needed anymore
    fullnames <- o <- NULL
    verbose && exit(verbose)
  }

  verbose <- less(verbose)
  names <- getNames(this)
  for (ii in arrays) {
    verbose && enter(verbose, sprintf("Array #%d of %d ('%s')",
                                         ii, nbrOfArrays, names[ii]))
    ce <- this[[ii]]

    verbose <- less(verbose, 50)
    verbose && enter(verbose, "Extracting estimates");  # 3-4s
    dataOne <- lapply(data, FUN=lapply, function(group) {
      # theta = group$theta[ii] = ...
      # stdvs = group$sdTheta[ii] = ...
      list(
        theta=.subset(.subset2(group, "theta"), ii),
        sdTheta=.subset(.subset2(group, "sdTheta"), ii),
        thetaOutliers=.subset(.subset2(group, "thetaOutliers"), ii)
      )
    })
    verbose && str(verbose, dataOne[1])
    verbose && exit(verbose)

    # Assert correct lengths  /HB 2007-04-12
    if (length(cdf) != length(dataOne)) {
      throw("Internal error: Lengths of 'cdf' and 'dataOne' differ: ",
                                    length(cdf), " != ", length(dataOne))
    }

    verbose && enter(verbose, "Updating file");  # 6-7s ~98% in encode()
#    verbose && printf(verbose, "class(ce)[1]: %s\n", class(ce)[1])
#    updateUnits(ce, cdf=cdf, data=dataOne, verbose=less(verbose, 50))
    updateUnits(ce, cdf=cdf, data=dataOne, verbose=verbose)
    # Not needed anymore
    dataOne <- ce <- NULL
    verbose && exit(verbose)

    verbose <- more(verbose, 50)

    gc <- gc()
    verbose && print(verbose, gc)

    verbose && exit(verbose)
  } # for (ii ...)
  verbose <- more(verbose)
  verbose && exit(verbose)
}, protected=TRUE)



setMethodS3("getAverageFile", "ChipEffectSet", function(this, indices="remaining", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'indices':
  if (identical(indices, "remaining")) {
  } else if (is.null(indices)) {
#    # Update only cells which stores values
#    indices <- getCellIndices(this, verbose=verbose)
#    indices <- unlist(indices, use.names=FALSE)
  }

  NextMethod("getAverageFile", indices=indices)
})


setMethodS3("findUnitsTodo", "ChipEffectSet", function(this, ...) {
  # Look into the chip-effect file that comes last in a lexicographic
  # order, becuase that is updated last.
  names <- getFullNames(this)
  idx <- order(names, decreasing=TRUE)[1]
  df <- this[[idx]]
  findUnitsTodo(df, ...)
})



setMethodS3("extractMatrix", "ChipEffectSet", function(this, ..., field=c("theta", "sdTheta", "RLE", "NUSE"), drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'field':
  field <- match.arg(field)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  if (toupper(field) == "RLE") {
    verbose && enter(verbose, "Extracting chip effects and calculating RLE scores")

    # RLE - Relative Log2 Expression values
    # Get chip effect estimates (on the log scale)
    data <- extractMatrix(this, ..., field="theta", verbose=less(verbose, 1))
    data <- log2(data);  # ...stored on the intensity scale

    # Robust average (on the log scale, but stored on the intensity scale!)
    avg <- getAverageLog(this, field="intensities", mean="median",
                                                      verbose=less(verbose,1))
    dataR <- extractMatrix(avg, field="theta", ..., verbose=less(verbose, 1))
    dataR <- log2(as.vector(dataR))

    # Log ratios of chip effects
    data <- data - dataR
    # Not needed anymore
    dataR <- avg <- NULL

    verbose && exit(verbose)
  } else if (toupper(field) == "NUSE") {
    verbose && enter(verbose, "Extracting standard errors for chip effects and calculates NUSE scores")
    # NUSE - Normalized Unscaled Standard Errors
    # Get standard errors (on the log scale)
    data <- extractMatrix(this, ..., field="sdTheta",
                                                   verbose=less(verbose, 1))
    data <- log2(data);  # ...stored on the intensity scale

    # Robust average (on the log scale, but stored on the intensity scale!)
    avg <- getAverageLog(this, field="stdvs", mean="median",
                                                    verbose=less(verbose, 1))
    dataR <- extractMatrix(avg, field="theta", ..., verbose=less(verbose, 1))
    dataR <- log2(as.vector(dataR))

    # Log ratios of standard errors
    data <- data / dataR
    # Not needed anymore
    dataR <- avg <- NULL

    verbose && exit(verbose)
  } else {
    data <- NextMethod("extractMatrix", field=field)
  }

  # Drop singleton dimensions?
  if (drop) {
    data <- drop(data)
  }

  data
})
