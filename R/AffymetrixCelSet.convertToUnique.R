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
  # getting output directory
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  inputTags <- paste( getTags(this), collapse=",")
  verbose && cat(verbose, "Input tags:", inputTags);
  
  curPath <- getPath(this)
  dataDir <- gsub("rawData","probeData",strsplit(curPath, "/")[[1]][1])
  allTags <- gsub("^,","",paste( inputTags, tags, sep="," ))
  chipType <- getChipType(this,fullname=FALSE)
  
  outputPath <- paste( dataDir, paste(getName(this),allTags,sep=","), chipType, sep="/" )
  verbose && cat(verbose, "Input Path: ", curPath);
  verbose && cat(verbose, "Output Path:", outputPath);
  verbose && cat(verbose, "allTags:", allTags);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # check if already done
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Test whether dataset exists");
  # HB: Don't think argument 'chipType' makes a difference if 'cdf' is given.
  outputDataSet <- NULL
  tryCatch({
    outputDataSet <- AffymetrixCelSet$byName(getName(this), tags=allTags, verbose=verbose, 
						  cdf=cdfUnique, paths=dataDir, checkChipType=FALSE);
  }, error = function(ex) {});
  
  if (inherits(outputDataSet, "AffymetrixCelSet")) {
    verbose && cat(verbose, "Dataset already created.");
    verbose && exit(verbose);
    return(invisible(outputDataSet));
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
      verbose && enter(verbose, paste("Converting CEL data from standard to unique CDF for sample ", kk, " ( ", getName(df), " ) of ", nbrOfArrays,sep=""));
  
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Read data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading intensity values according to standard CDF");
      data <- readCelUnits(getPathname(df), cdf=cdfStandard, dropArrayDim=TRUE);
      #hdr <- readCelHeader(getPathname(df));
      verbose && exit(verbose);

      fullname <- getFullName(df);
      filename <- sprintf("%s.CEL", fullname);
      pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);
	  
	  # Build a valid CEL header
      celHeader <- cdfHeaderToCelHeader(cdfHeader, sampleName=fullname);

      # Add some extra information about what the CEL file is for
      params <- c(Descripion="This CEL file was created by the aroma.affymetrix package.");
      parameters <- gsub(" ", "_", params);
      names(parameters) <- names(params);
      parameters <- paste(names(parameters), parameters, sep=":");
      parameters <- paste(parameters, collapse=";");
      parameters <- paste(celHeader$parameters, parameters, "", sep=";");
      parameters <- gsub(";;", ";", parameters);
      parameters <- gsub(";$", "", parameters);
      celHeader$parameters <- parameters;

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Write data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");
      createCel(pathname, header=celHeader);
      verbose && cat(verbose, "Writing values according to unique CDF");
      updateCelUnits(pathname, cdf=cdfUniqueIndices, data=data, verbose=FALSE);
      verbose && exit(verbose);

      rm(data);
      gc <- gc();
      verbose && print(verbose, gc);

      verbose && exit(verbose);
  } # for (kk ...)

  outputDataSet <- AffymetrixCelSet$byName(getName(this), tags=allTags, 
                      cdf=cdfUnique, verbose=verbose, checkChipType=FALSE);
  verbose && exit(verbose);
  
  invisible(outputDataSet);
})


############################################################################
# HISTORY:
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
