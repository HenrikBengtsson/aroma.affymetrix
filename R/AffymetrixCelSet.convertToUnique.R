setMethodS3("convertToUnique", "AffymetrixCelSet", function(this, ..., tags="UNQ", force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Converting to unique CDF");
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already unique?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this)
  if (isUniqueCdf(cdf)) {
    verbose && cat(verbose, "Already based on a unique CDF");
    verbose && exit(verbose);
    return(invisible(this));
  } else {
    verbose && enter(verbose, "Getting unique CDF");
    cdfUnique <- getUniqueCdf(cdf)
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting output directory
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rootPath <- "probeData";
  verbose && cat(verbose, "Output root: ", rootPath);

  srcTags <- getTags(this, collapse=",");
  verbose && cat(verbose, "Source tags: ", srcTags);

  verbose && cat(verbose, "User tags: ", tags);
  
  tags <- c(srcTags, tags);
  tags <- tags[nchar(tags) > 0];
  tags <- paste(tags, collapse=",");
  verbose && cat(verbose, "Output tags: ", tags);

  chipType <- getChipType(this, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);
  
  fullname <- paste(c(getName(this), tags), collapse=",");
  verbose && cat(verbose, "Output fullname: ", fullname);

  outputPath <- file.path(rootPath, fullname, chipType);
  outputPath <- Arguments$getWritablePath(outputPath);
  verbose && cat(verbose, "Output Path: ", outputPath);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check if already done
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Test whether dataset exists");
  # HB: Don't think argument 'chipType' makes a difference if 'cdf' is given.
  outputDataSet <- NULL
  tryCatch({
    res <- AffymetrixCelSet$byName(fullname, cdf=cdfUnique, 
                     paths=rootPath, checkChipType=FALSE, verbose=verbose);
  }, error = function(ex) {});
  
  if (inherits(res, "AffymetrixCelSet")) {
    srcFullnames <- getFullNames(this);
    fullnames <- getFullNames(res);
    missing <- setdiff(srcFullnames, fullnames);
    if (length(missing) == 0) {
      verbose && cat(verbose, "Dataset already created.");
      verbose && exit(verbose);
      return(invisible(res));
    }
  }
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read indices for old and new
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading cell indices from standard CDF");
  cdfStandard <- readCdf(getPathname(cdf), units=NULL, readXY=FALSE, readBases=FALSE, readIndexpos=FALSE, readAtoms=FALSE,readUnitType=FALSE, readUnitDirection=FALSE, readUnitNumber=FALSE, readUnitAtomNumbers=FALSE, readGroupAtomNumbers=FALSE, readGroupDirection=FALSE, readIndices=TRUE, readIsPm=FALSE);
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Reading cell indices list from unique CDF");
  cdfUniqueIndices <- readCdf(getPathname(cdfUnique), units=NULL, readXY=FALSE, readBases=FALSE, readIndexpos=FALSE, readAtoms=FALSE,readUnitType=FALSE, readUnitDirection=FALSE, readUnitNumber=FALSE, readUnitAtomNumbers=FALSE, readGroupAtomNumbers=FALSE, readGroupDirection=FALSE, readIndices=TRUE, readIsPm=FALSE);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize all arrays simultaneously
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- nbrOfArrays(this);
  
  # Get CDF header
  cdfHeader <- getHeader(cdfUnique);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Do the conversion from standard CDF to unique CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq_len(nbrOfArrays)) {
      df <- getFile(this, kk);
      verbose && enter(verbose, sprintf("Converting CEL data from standard to unique CDF for sample #%d (%s) of %d", kk, getName(df), nbrOfArrays));
  
      dfFullname <- getFullName(df);
      filename <- sprintf("%s.CEL", dfFullname);
      pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Read data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading intensity values according to standard CDF");
      data <- readCelUnits(getPathname(df), cdf=cdfStandard, dropArrayDim=TRUE);
      verbose && exit(verbose);

      # Build a valid CEL header
      celHeader <- cdfHeaderToCelHeader(cdfHeader, sampleName=dfFullname);
      rm(dfFullname);

      # Add some extra information about what the CEL file is for
      params <- c(Descripion="This CEL file was created by the aroma.affymetrix package.");
      parameters <- gsub(" ", "_", params, fixed=TRUE);
      names(parameters) <- names(params);
      parameters <- paste(names(parameters), parameters, sep=":");
      parameters <- paste(parameters, collapse=";");
      parameters <- paste(celHeader$parameters, parameters, "", sep=";");
      parameters <- gsub(";;", ";", parameters, fixed=TRUE);
      parameters <- gsub(";$", "", parameters);
      celHeader$parameters <- parameters;

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Write data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results");

      # Write to a temporary file
      pathnameT <- sprintf("%s.tmp", pathname);
      verbose && cat(verbose, "Temporary pathname: ", pathnameT);

      createCel(pathnameT, header=celHeader);
      verbose && cat(verbose, "Writing values according to unique CDF");
      updateCelUnits(pathnameT, cdf=cdfUniqueIndices, data=data, verbose=FALSE);
      verbose && exit(verbose);

      rm(data);
      gc <- gc();
      verbose && print(verbose, gc);

      # Rename temporary file
      verbose && enter(verbose, "Renaming temporary file");
      res <- file.rename(pathnameT, pathname);
      if (!isFile(pathname)) {
        throw("Failed to rename temporary file (final file does not exist): ", pathnameT, " -> ", pathname);
      }
      if (isFile(pathnameT)) {
        throw("Failed to rename temporary file (temporary file still exists): ", pathnameT, " -> ", pathname);
      }
      rm(pathnameT);
      verbose && exit(verbose);

      verbose && exit(verbose);
  } # for (kk ...)

  res <- AffymetrixCelSet$byName(fullname, cdf=cdfUnique,
                                      checkChipType=FALSE, verbose=verbose);
  verbose && exit(verbose);
  
  invisible(res);
})


############################################################################
# HISTORY:
# 2010-05-13 [HB]
# o Yday's fixes had some hiccups.
# 2010-05-12 [HB]
# o BUG FIX: convertToUnique() for AffymetrixCelSet would not recognize
#   Windows Shortcut links.
# o ROBUSTNESS: Now convertToUnique() for AffymetrixCelSet writes to 
#   temporary files which are renamed when complete.  This lowers the
#   risk of generating incomplete files.
# o CLEAN UP: Code cleanup.
# 2009-03-18 [MR]
# o changed the way CEL headers are made ... it now uses cdfHeaderToCelHeader()
# 2008-12-08 [MR]
# o fixed small bug when operating on raw data
# 2008-12-04 [MR]
# o fixed small bug when previous dataset not available
# 2008-11-28 [HB]
# o Tidying up code.
# o Replaced try() with tryCatch() statement.  I consider try() to be 
#   obsolete.
# 2008-11-11 [MR]
# o Created.
############################################################################ 
