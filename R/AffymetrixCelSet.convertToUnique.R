# @author "MR, HB"
setMethodS3("convertToUnique", "AffymetrixCelSet", function(this, ..., tags="*", force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (is.null(tags)) tags <- c("*")
  tags <- Arguments$getTags(tags, collapse=NULL)
  tags[tags == "*"] <- "UNQ"

  # Argument 'force':
  force <- Arguments$getLogical(force)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Converting to unique CDF")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already unique?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this)

  # Already done?
  if (isUniqueCdf(cdf)) {
    verbose && cat(verbose, "Already based on a unique CDF")
    verbose && exit(verbose)
    return(invisible(this))
  }

  verbose && enter(verbose, "Getting unique CDF")
  cdfUnique <- getUniqueCdf(cdf)
  verbose && exit(verbose)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting output directory
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rootPath <- "probeData"
  verbose && cat(verbose, "Output root: ", rootPath)

  srcTags <- getTags(this, collapse=",")
  verbose && cat(verbose, "Source tags: ", srcTags)

  verbose && cat(verbose, "User tags: ", tags)

  tags <- c(srcTags, tags)
  tags <- tags[nzchar(tags)]
  tags <- paste(tags, collapse=",")
  verbose && cat(verbose, "Output tags: ", tags)

  chipType <- getChipType(this, fullname=FALSE)
  verbose && cat(verbose, "Chip type: ", chipType)

  fullname <- paste(c(getName(this), tags), collapse=",")
  verbose && cat(verbose, "Output fullname: ", fullname)

  outputPath <- file.path(rootPath, fullname, chipType)
  outputPath <- Arguments$getWritablePath(outputPath)
  verbose && cat(verbose, "Output Path: ", outputPath)


  # Expected output set
  fullnames <- getFullNames(this)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Checking if dataset already exists")

  res <- tryCatch({
    # HB: Don't think argument 'checkChipType' makes a difference if
    #     argument 'cdf' is given.
    AffymetrixCelSet$byName(fullname, cdf=cdfUnique,
                           checkChipType=FALSE, verbose=less(verbose, 10))
  }, error = function(ex) { NULL })

  if (inherits(res, "AffymetrixCelSet")) {
    # Extract samples in the same order as they appear in the input
    # data set, and if more were found, drop those.
    # WORKAROUND/TO BE REMOVED: R.utils (<= 1.19.0) will give an error
    # with extract() if length(fullnames) > length(res). /HB 2012-12-01
    if (length(res) >= length(fullnames)) {
      res <- extract(res, fullnames, onMissing="drop", onDuplicates="error")
    }

    # Is output set complete?
    missing <- setdiff(fullnames, getFullNames(res))
    if (length(missing) == 0) {
      if (!force) {
        verbose && cat(verbose, "Detected existing output dataset. Skipping.")
        verbose && exit(verbose)
        verbose && exit(verbose)
        return(invisible(res))
      }
      verbose && cat(verbose, "Detected existing output dataset, but will force reprocessing.")
    } else if (length(missing) > 0) {
      verbose && cat(verbose, "Detected partial output dataset.")
    }
  }

  # Not needed anymore
  res <- NULL;   # Not needed anymore
  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read indices for old and new
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading cell indices from standard CDF")
  cdfStandard <- .readCdf(getPathname(cdf), units=NULL, readXY=FALSE, readBases=FALSE, readIndexpos=FALSE, readAtoms=FALSE,readUnitType=FALSE, readUnitDirection=FALSE, readUnitNumber=FALSE, readUnitAtomNumbers=FALSE, readGroupAtomNumbers=FALSE, readGroupDirection=FALSE, readIndices=TRUE, readIsPm=FALSE)
  verbose && exit(verbose)

  verbose && enter(verbose, "Reading cell indices list from unique CDF")
  cdfUniqueIndices <- .readCdf(getPathname(cdfUnique), units=NULL, readXY=FALSE, readBases=FALSE, readIndexpos=FALSE, readAtoms=FALSE,readUnitType=FALSE, readUnitDirection=FALSE, readUnitNumber=FALSE, readUnitAtomNumbers=FALSE, readGroupAtomNumbers=FALSE, readGroupDirection=FALSE, readIndices=TRUE, readIsPm=FALSE)

  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Process all arrays simultaneously
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(this)

  # Get CDF header
  cdfHeader <- getHeader(cdfUnique)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Do the conversion from standard CDF to unique CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- listenv()

  for (kk in seq_len(nbrOfArrays)) {
      df <- this[[kk]]
      verbose && enter(verbose, sprintf("Converting CEL data from standard to unique CDF for sample #%d (%s) of %d", kk, getName(df), nbrOfArrays))

      dfFullname <- getFullName(df)
      filename <- sprintf("%s.CEL", dfFullname)
      pathname <- Arguments$getWritablePathname(filename, path=outputPath, mustNotExist=FALSE, ...)

      # Skip?
      isFile <- isFile(pathname)
      if (!force && isFile) {
        verbose && cat(verbose, "Already processed. Skipping.")
        res[[kk]] <- pathname
        verbose && exit(verbose)
        next
      }


      res[[kk]] %<-% {
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Read data
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        verbose && enter(verbose, "Reading intensity values according to standard CDF")
        data <- .readCelUnits(getPathname(df), cdf=cdfStandard, dropArrayDim=TRUE)
        verbose && exit(verbose)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Write data
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        verbose && enter(verbose, "Creating CEL header")
        # Build a valid CEL header
        celHeader <- .cdfHeaderToCelHeader(cdfHeader, sampleName=dfFullname)

        # Add some extra information about what the CEL file is for
        params <- c(Descripion="This CEL file was created by the aroma.affymetrix package.")
        parameters <- gsub(" ", "_", params, fixed=TRUE)
        names(parameters) <- names(params)
        parameters <- paste(names(parameters), parameters, sep=":")
        parameters <- paste(parameters, collapse=";")
        parameters <- paste(celHeader$parameters, parameters, "", sep=";")
        parameters <- gsub(";;", ";", parameters, fixed=TRUE)
        parameters <- gsub(";$", "", parameters)
        celHeader$parameters <- parameters
        verbose && exit(verbose)

        # Create CEL file to store results, if missing
        verbose && enter(verbose, "Creating CEL file for results")

        # Remove existing file (=overwrite)?
        if (isFile) file.remove(pathname)

        # Write to a temporary file
        pathnameT <- pushTemporaryFile(pathname, verbose=verbose)

        .createCel(pathnameT, header=celHeader)
        verbose && cat(verbose, "Writing values according to unique CDF")
        .updateCelUnits(pathnameT, cdf=cdfUniqueIndices, data=data, verbose=FALSE)
        verbose && exit(verbose)

        # Not needed anymore
        gc <- gc()
        verbose && print(verbose, gc)

        # Rename temporary file
        popTemporaryFile(pathnameT, verbose=verbose)

        ## Create checksum file
        dfZ <- getChecksumFile(pathname)

        pathname
      } ## %<-%

      verbose && exit(verbose)
  } # for (kk ...)

  ## Resolve futures
  res <- as.list(res)
  res <- NULL

  res <- AffymetrixCelSet$byName(fullname, cdf=cdfUnique,
                                      checkChipType=FALSE, verbose=verbose)

  # Extract samples in the same order as they appear in the input
  # data set, and if more were found, drop those.
  res <- extract(res, fullnames, onMissing="error", onDuplicates="error")

  verbose && exit(verbose)

  invisible(res)
}, protected=TRUE) # convertToUnique()
