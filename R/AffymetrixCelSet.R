###########################################################################/**
# @RdocClass AffymetrixCelSet
#
# @title "The AffymetrixCelSet class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixCelSet object represents a set of Affymetrix CEL files 
#  with \emph{identical} chip types.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{files}{A @list of @see "AffymetrixCelFile":s.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \examples{\dontrun{
#   @include "../incl/AffymetrixCelSet.Rex"
# }}
#
# \seealso{
#   @see "AffymetrixCelFile".
# }
#
# @author
#*/###########################################################################
setConstructorS3("AffymetrixCelSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    reqFileClass <- "AffymetrixCelFile";
    base::lapply(files, FUN=function(df) {
      if (!inherits(df, reqFileClass))
        throw("Argument 'files' contains a non-", reqFileClass, 
                                                  " object: ", class(df)[1]);
    })
  } else if (inherits(files, "AffymetrixCelSet")) {
    return(as.AffymetrixCelSet(files));
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }


  this <- extend(AffymetrixFileSet(files=files, ...), "AffymetrixCelSet",
    "cached:.intensities" = NULL,
    "cached:.intensitiesIdxs" = NULL,
    "cached:.readUnitsCache" = NULL,
    "cached:.getUnitIntensitiesCache" = NULL,
    "cached:.averageFiles" = list(),
    "cached:.timestamps" = NULL,
    "cached:.fileSize" = NULL
  );

  if (length(this$files) > 0) {
    # Make sure the set name is non-empty
    name <- getName(this);
    if (nchar(name) == 0) {
      throw("An ", class(this)[1], " must have a name of at least length one: ", this$.pathname);
    }
  }

  this;
})


setMethodS3("clearCache", "AffymetrixCelSet", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".intensities", ".intensitiesIdxs", ".readUnitsCache", 
           ".getUnitIntensitiesCache", ".timestamps", ".fileSize")) {
    this[[ff]] <- NULL;
  }
  this$.averageFiles <- list();

  if (nbrOfFiles(this) > 0) {
    # Clear the cache for the CDF.
    cdf <- getCdf(this);
    clearCache(cdf);
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)


setMethodS3("clone", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Cloning Affymetrix CEL set");

  # Clone itself and the files.  The call below will clear the cache!
  object <- NextMethod("clone", clear=TRUE, ..., verbose=less(verbose));
  clearCache(object);

  if (nbrOfFiles(object) > 0) {
    # Clone the CDF (this will update the CDF of all file object)
    verbose && enter(verbose, "Cloning CDF");
    cdf <- getCdf(object);
    cdf <- clone(cdf);
    verbose && exit(verbose);
    verbose && enter(verbose, "Adding CDF to CEL set");
    setCdf(object, cdf, .checkArgs=FALSE);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  object;
}, private=TRUE)



setMethodS3("append", "AffymetrixCelSet", function(this, other, clone=TRUE, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (!inherits(other, class(this)[1])) {
    throw("Argument 'other' is not an ", class(this)[1], " object: ", 
                                                      class(other)[1]);
  }

  verbose && enter(verbose, "Appending CEL set");
  verbose && print(verbose, other);

  # Validate chip type
  cdf <- getCdf(this);
  chipType <- getChipType(cdf);
  for (file in getFiles(other)) {
    oCdf <- getCdf(file);
    oChipType <- getChipType(oCdf);
    if (!identical(oChipType, chipType)) {
      throw("Argument 'other' contains a CEL file of different chip type: ",
                                                oChipType, " != ", chipType);
    }
  }

  # Append other
  this <- NextMethod("append", this, other=other, clone=clone, ...);

  # Set the same CDF for all CEL files
  verbose && enter(verbose, "Updating the CDF for all files");
  setCdf(this, cdf);
  verbose && exit(verbose);

  verbose && exit(verbose);

  this;
})




###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the Affymetrix CEL set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("as.character", "AffymetrixCelSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  s <- c(s, sprintf("Tags: %s", tags));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Platform: %s", getPlatform(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  n <- nbrOfArrays(this);
  s <- c(s, sprintf("Number of arrays: %d", n));
  names <- getNames(this);
  if (n >= 5)
    names <- c(names[1:2], "...", names[n]);
  names <- paste(names, collapse=", ");
  s <- c(s, sprintf("Names: %s", names));

  # Get CEL header timestamps?
  maxCount <- getOption(aromaSettings, "output/timestampsThreshold");
  if (maxCount >= nbrOfArrays(this)) {
    ts <- getTimestamps(this);
    # Note: If ts <- range(ts) is used and the different timestamps uses
    # tifferent 'tzone' attributes, e.g. if some files where scanning during
    # daylight savings time and some not, we will get a warning saying:
    # "'tzone' attributes are inconsistent".  By doing the below, we avoid
    # this warning (which confuses users). 
    ts <- sort(ts);
    ts <- ts[c(1,n)];
    ts <- format(ts, "%Y-%m-%d %H:%M:%S");  # range() gives strange values?!?
    s <- c(s, sprintf("Time period: %s -- %s", ts[1], ts[2]));
  } else {
    s <- c(s, sprintf("Time period: [not reported if more than %.0f arrays]", as.double(maxCount)));
  }
  s <- c(s, sprintf("Total file size: %.2fMB", getFileSize(this)/1024^2));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getTimestamps", "AffymetrixCelSet", function(this, ..., force=FALSE) {
  ts <- this$.timestamps;

  if (force || is.null(ts)) {
    # Get CEL header dates
    ts <- lapply(this, getTimestamp);
    ts <- do.call("c", args=ts);
    this$.timestamps <- ts;
  }

  ts;
})


setMethodS3("getIdentifier", "AffymetrixCelSet", function(this, ..., force=FALSE) {
  identifier <- this$.identifier;
  if (force || is.null(identifier)) {
    identifier <- NextMethod("getIdentifier");
    if (is.null(identifier)) {
      identifiers <- lapply(this, getIdentifier);
      identifier <- digest2(identifiers);
    }
    this$.identifier <- identifier;
  }
  identifier;
}, private=TRUE)



###########################################################################/**
# @RdocMethod getChipType
#
# @title "Gets the chip type for this CEL set"
#
# \description{
#  @get "title".
#  The chip type is inferred from the current CDF.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @character string.
# }
#
# @author
#
# \seealso{
#   @seemethod "getCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getChipType", "AffymetrixCelSet", function(this, ...) {
  unf <- getUnitNamesFile(this);
  getChipType(unf, ...);
})

setMethodS3("getPlatform", "AffymetrixCelSet", function(this, ...) {
  "Affymetrix";
})


setMethodS3("getUnitNamesFile", "AffymetrixCelSet", function(this, ...) {
  aFile <- getFile(this, 1);
  getUnitNamesFile(aFile, ...);
})



###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this CEL set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCdfFile" object.
# }
#
# @author
#
# \seealso{
#   @seemethod "setCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCdf", "AffymetrixCelSet", function(this, ...) {
  getCdf(this$files[[1]], ...);
})


###########################################################################/**
# @RdocMethod setCdf
#
# @title "Sets the CDF structure for this CEL set"
#
# \description{
#  @get "title".  This structures is assigned to all CEL files in the set.
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{An @see "AffymetrixCdfFile" object.}
#   \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#     May also be a @see "R.utils::Verbose" object.}
#   \item{...}{Not used.}
#   \item{.checkArgs}{(Internal) If @FALSE, arguments are not validated.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "getCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("setCdf", "AffymetrixCelSet", function(this, cdf, verbose=FALSE, ..., .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    # Argument 'cdf':
    if (!inherits(cdf, "AffymetrixCdfFile")) {
      throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
    }
  
    # Assure that the CDF is compatible with the CEL file
    if (nbrOfFiles(this) > 0) {
      cf <- getFile(this, 1);
      if (nbrOfCells(cdf) != nbrOfCells(cf)) {
        throw("The specified CDF structure ('", getChipType(cdf), "') is not compatible with the chip type ('", getChipType(cf), "') of the CEL file. The number of cells do not match: ", nbrOfCells(cdf), " != ", nbrOfCells(cf));
      }
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Setting CDF for CEL set");
  verbose && print(verbose, cdf);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # ASCII CDFs are only allowed if explicited accepted in the rules.
  # To change, do:
  #  setSetting(aroma.affymetrix, "rules$allowAsciiCdfs", TRUE)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  allowAsciiCdfs <- getOption(aromaSettings, "rules/allowAsciiCdfs", FALSE);
  if (allowAsciiCdfs) {
    # ASCII CDF are *not* allowed
    ff <- getFileFormat(cdf);
    if (regexpr("ASCII", ff) != -1) {
      throw("Cannot set CDF for data set. The given CDF is in ASCII format, which is protected against by default. It is much faster and more memory efficient to work with binary CDF files. Use affxparser::convertCdf() to convert a CDF into another format.  Use setSetting(aroma.affymetrix, \"rules$allowAsciiCdfs\", TRUE) to allow ASCII CDFs. For more details, see the online help pages. Details on the CDF file: ", getPathname(cdf), " [", ff, "].");
    }
  }

  # Set the CDF for all CEL files
  verbose && enter(verbose, "Setting CDF for each CEL file");
  lapply(this, setCdf, cdf, .checkArgs=FALSE, ...);
  verbose && exit(verbose);

  # Have to clear the cache 
  verbose && enter(verbose, "Clearing data-set cache");
  clearCache(this);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(this);
})


setMethodS3("findByName", "AffymetrixCelSet", function(static, name, tags=NULL, chipType=NULL, paths=c("rawData", "probeData"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'name':
  name <- Arguments$getCharacter(name, length=c(1,1));
  if (nchar(name) == 0) {
    throw("A ", class(static)[1], " must have a non-empty name: ''");
  }

  # Arguments 'paths':
  if (is.null(paths)) {
    paths <- eval(formals(findByName.AffymetrixCelSet)[["paths"]]);
  }


  # Look only in existing directories
  paths <- sapply(paths, FUN=filePath, expandLinks="any");
  paths0 <- paths;
  paths <- paths[sapply(paths, FUN=isDirectory)];
  if (length(paths) == 0) {
    throw("None of the data directories exist: ", 
                                           paste(paths0, collapse=", "));
  }

  # The full name of the data set
  fullname <- paste(c(name, tags), collapse=",");

  # Look for matching data sets
  paths <- file.path(paths, fullname);
  paths <- paths[sapply(paths, FUN=isDirectory)];
  if (length(paths) == 0)
    return(NULL);

  # Look for matching chip type sets?
  if (!is.null(chipType)) {
    paths <- file.path(paths, chipType);
    paths <- paths[sapply(paths, FUN=isDirectory)];
    if (length(paths) == 0)
      return(NULL);
  }

  if (length(paths) > 1) {
    warning("Found duplicated data set: ", paste(paths, collapse=", "));
    paths <- paths[1];
  }
  
  paths;
}, static=TRUE)


setMethodS3("fromName", "AffymetrixCelSet", function(static, ...) {
  byName(static, ...);
}, static=TRUE)

setMethodS3("byName", "AffymetrixCelSet", function(static, name, tags=NULL, chipType=NULL, cdf=NULL, paths=NULL, ...) {
  # Argument 'cdf':
  if (!is.null(cdf)) {
    if (!inherits(cdf, "AffymetrixCdfFile"))
      throw("Argument 'cdf' must be an AffymetrixCdfFile object: ", class(cdf)[1]);
  }

  # Argument 'chipType':
  if (is.null(chipType)) {
    if (!is.null(cdf)) {
      chipType <- getChipType(cdf, fullname=FALSE);  # Without tags
    } else {
      throw("Argument 'chipType' must be specified unless argument 'cdf' is specified.");
    }
  }

  suppressWarnings({
    path <- findByName(static, name, tags=tags, chipType=chipType, paths=paths, ...);
  })
  if (is.null(path)) {
    path <- file.path(paste(c(name, tags), collapse=","), chipType);
    throw("Cannot create ", class(static)[1], ".  No such directory: ", path);
  }

  suppressWarnings({
    fromFiles(static, path=path, cdf=cdf, ...);
  })
}, static=TRUE)


setMethodS3("update2", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating ", class(this)[1]);
  updateSampleAnnotationSet(this, ..., verbose=less(verbose, 1));
  verbose && exit(verbose);
}, protected=TRUE)


setMethodS3("updateSampleAnnotationSet", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Scan for SAF files and apply them
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Scanning for and applying sample annotation files");
  sasPath <- "annotationData/samples/";
  sasPath <- filePath(sasPath, expandLinks="any");
  mkdirs(sasPath);

  sas <- SampleAnnotationSet$fromPath(sasPath, verbose=less(verbose));
  if (nbrOfFiles(sas) == 0) {
    verbose && cat(verbose, "No sample annotation files found.");
  } else {
    verbose && print(verbose, sas);
    setAttributesBy(this, sas);
  }
  # Store the SAFs for now.
  this$.sas <- sas;

  verbose && exit(verbose);

  invisible(this);
}, private=TRUE)  # updateSampleAnnotationSet()


setMethodS3("fromFiles", "AffymetrixCelSet", function(static, path="rawData/", pattern="[.](c|C)(e|E)(l|L)$", cdf=NULL, checkChipType=is.null(cdf), ..., onDuplicates=c("keep", "exclude", "error"), fileClass="AffymetrixCelFile", force=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'onDuplicates':
  onDuplicates <- match.arg(onDuplicates);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  
  verbose && enter(verbose, "Defining ", class(static)[1], " from files");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Look for cached results (useful for extremely large data set)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="fromFiles", class=class(static)[1], path=path, pattern=pattern, cdf=cdf, checkChipType=checkChipType, ..., fileClass=fileClass);
  dirs <- "aroma.affymetrix";
  res <- loadCache(key=key, dirs=dirs);
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "Found cached results");
    verbose && exit(verbose);
    return(res);
  }

  set <- fromFiles.AffymetrixFileSet(static, path=path, pattern=pattern, ..., fileClass=fileClass, verbose=less(verbose));

  verbose && enter(verbose, "Retrieved files: ", nbrOfFiles(set));

  if (nbrOfFiles(set) > 0) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Handle duplicates
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if (onDuplicates %in% c("exclude", "error")) {
      dups <- isDuplicated(set, verbose=less(verbose));
      ndups <- sum(dups);
      if (ndups > 0) {
        dupsStr <- paste(getNames(set)[dups], collapse=", ");
        if (onDuplicates == "error") {
          msg <- paste("Detected ", ndups, " duplicated CEL files (same datestamp): ", dupsStr, sep="");
          verbose && cat(verbose, "ERROR: ", msg);
          throw(msg);
        } else if (onDuplicates == "exclude") {
          set <- extract(set, !dups);
          msg <- paste("Excluding ", ndups, " duplicated CEL files (same datestamp): ", dupsStr, sep="");
          verbose && cat(verbose, "WARNING: ", msg);
          warning(msg);
        }
      }
    }
  
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Scan all CEL files for possible chip types
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Chip type according to the directory structure
    path <- getPath(set);
    chipType <- basename(path);
    verbose && cat(verbose, 
                   "The chip type according to the path is: ", chipType);
  
    # Let the directory name specify the chip type?
    if (checkChipType) {
      verbose && enter(verbose, "Scanning CEL set for used chip types");
      # This takes time if the CEL files are ASCII files.
      chipTypes <- sapply(set, FUN=function(file) {
        readCelHeader(getPathname(file))$chiptype;
      })
      tChipTypes <- table(chipTypes);
      verbose && print(verbose, tChipTypes);
      nbrOfChipTypes <- length(tChipTypes);
      verbose && exit(verbose);
    
      if (nbrOfChipTypes > 1) {
        verbose && cat(verbose, "Detected ", nbrOfChipTypes, 
                                        " different chip types in CEL set: ",
                                    paste(names(tChipTypes), collapse=", "));
      } else {
        verbose && cat(verbose, "All CEL files use the same chip type: ", 
                                                          names(tChipTypes));
      }
    
      # If chip type is taken from CEL headers and there are more than
      # one chip type in the set, then it is an error.
      if (nbrOfChipTypes > 1) {
        throw("Detected ", nbrOfChipTypes, " different chip types in CEL set. Use argument 'checkChipType=FALSE' to let the directory name of the CEL set specify the chip type instead: ", paste(names(tChipTypes), collapse=", "));
      } 
  
      # Validate that the directory name matches the chip type
      tChipTypesShort <- gsub(",.*", "", names(tChipTypes));
      if (!identical(tChipTypesShort, chipType)) {
        throw("Invalid name of directory containing CEL files. The name of the directory (", chipType, ") must be the same as the chip type used for the CEL files (", names(tChipTypes), ") unless using argument 'checkChipType=FALSE': ", path);
      }
    } else {
      verbose && cat(verbose, "Since 'checkChipType=FALSE', then the chip type specified by the directory name is used: ", chipType);
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Use the same CDF object for all CEL files.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if (is.null(cdf)) {
      verbose && enter(verbose, "Retrieving the CDF for chip type '", chipType, "' inferred from path");
      cf <- getFile(set, 1);
      nbrOfCells <- nbrOfCells(cf);
      cdf <- AffymetrixCdfFile$byChipType(chipType, nbrOfCells=nbrOfCells);
      verbose && exit(verbose);
  
      verbose && enter(verbose, "Check compatibility with 1st CEL file");
      verbose && cat(verbose, "Chip type: ", chipType);
      if (nbrOfCells(cdf) != nbrOfCells(cf)) {
        cdf <- getCdf(cf);
        chipType <- getChipType(cdf);
        verbose && cat(verbose, "Chip type (updated): ", chipType);
      }
      verbose && exit(verbose);
    } else {
      verbose && cat(verbose, "Using prespecified CDF: ", 
                     getChipType(cdf, fullname=TRUE));
    }
  }

  verbose && enter(verbose, "Updating the CDF for all files");
  setCdf(set, cdf);
  verbose && exit(verbose);

  # Let the new CEL set update itself
  update2(set, verbose=less(verbose, 1));

  # Save to file cache
  saveCache(set, key=key, dirs=dirs);

  verbose && exit(verbose);

  set;
}, protected=TRUE, static=TRUE)




###########################################################################/**
# @RdocMethod nbrOfArrays
#
# @title "Gets the number of arrays in the file set"
#
# \description{
#   @get "title".
#   This is just a wrapper for \code{nbrOfFiles()}.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfArrays", "AffymetrixCelSet", function(this, ...) {
  nbrOfFiles(this, ...);
})



###########################################################################/**
# @RdocMethod as.AffymetrixCelSet
# @alias as.AffymetrixCelSet.list
# @alias as.AffymetrixCelSet.default
#
# @title "Coerce an object to an AffymetrixCelSet object"
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
#   Returns an @see "AffymetrixCelSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.AffymetrixCelSet", "AffymetrixCelSet", function(object, ...) {
  object;
})

setMethodS3("as.AffymetrixCelSet", "list", function(object, ...) {
  AffymetrixCelSet(object, ...);
})

setMethodS3("as.AffymetrixCelSet", "default", function(object, ...) {
  throw("Cannot coerce object to an AffymetrixCelSet object: ", mode(object));
})



###########################################################################/**
# @RdocMethod isDuplicated
#
# @title "Identifies duplicated CEL files"
#
# \description{
#   @get "title" by comparing the timestamps in the CEL headers.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a @logical @vector of length equal to the number of files
#   in the set.
#   An element with value @TRUE indicates that the corresponding CEL file
#   has the same time stamp as another preceeding CEL file.
# }
#
# \examples{\dontrun{
#   # The data set of interest
#   ds <- AffymetrixCelSet$fromFiles(path=...)
#
#   # Added other data sets to be used as a reference
#   for (path in refPaths) {
#     dsR <- AffymetrixCelSet$fromFiles(path=path)
#     append(ds, dsR)
#   }
#
#   # Keep only unique arrays
#   ds <- extract(ds, !isDuplicated(ds))
# }}
#
# @author
#
# \seealso{
#   Internally @see "base::duplicated" is used to compare timestamps.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isDuplicated", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Scanning for duplicated data files");

  verbose && enter(verbose, "Reading the timestamp of all files");
  # Get the CEL header timestamp for all files
  timestamps <- getTimestamps(this);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying duplicated timestamps");
  dups <- duplicated(timestamps);
  names(dups) <- getNames(this);
  verbose && exit(verbose);

  if (verbose) {
    if (any(dups)) {
      cat(verbose, "Duplicated files:", paste(names(dups), collapse=", "));
    } else {
      cat(verbose, "Duplicated files: none.");
    }
  }

  verbose && exit(verbose);

  dups;
})



setMethodS3("getData", "AffymetrixCelSet", function(this, indices=NULL, fields=c("x", "y", "intensities", "stdvs", "pixels"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'indices':
  nbrOfCells <- nbrOfCells(getCdf(this));
  if (!is.null(indices)) {
    indices <- Arguments$getIndices(indices, range=c(1, nbrOfCells));
    nbrOfCells <- length(indices);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  nbrOfArrays <- nbrOfArrays(this);
  verbose && enter(verbose, "Getting cell data for ", nbrOfArrays, " arrays.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocating the return structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocating the return structure");
  nbrOfFields <- length(fields);
  res <- vector("list", nbrOfFields);
  names(res) <- fields;
  for (field in fields) {
    if (field %in% c("x", "y", "pixels")) {
      naValue <- as.integer(NA);
    } else {
      naValue <- as.double(NA);
    }
    res[[field]] <- matrix(naValue, nrow=nbrOfCells, ncol=nbrOfArrays);
  }
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading cell signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving data from ", nbrOfArrays, " arrays");
  for (kk in seq(length=nbrOfArrays)) {
    verbose && enter(verbose, "Array #", kk, " of ", nbrOfArrays);
    dataFile <- this$files[[kk]];
    value <- getData(dataFile, indices=indices, fields=fields, verbose=less(verbose));
    for (field in fields) {
      res[[field]][,kk] <- value[[field]];
      value[[field]] <- NULL;
    }
    rm(value); gc();
    verbose && exit(verbose);
  }
  verbose && exit(verbose);


  verbose && exit(verbose);

  res;
}) # getData()


###########################################################################/**
# @RdocMethod getIntensities
#
# @title "Gets cell intensities from a set of cells and a set of arrays"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Passed to @seemethod "getData".}
# }
#
# \value{
#   Returns a @numeric \eqn{NxK} matrix, where \eqn{N} is the number of 
#   cells read, and \eqn{K} is the number of arrays in the data set.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getIntensities", "AffymetrixCelSet", function(this, ...) {
  getData(this, ..., fields="intensities")$intensities;
}) # getIntensities()



###########################################################################/**
# @RdocMethod getUnitIntensities
#
# @title "Gets cell signals for a subset of units and a subset of arrays"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{units}{An @integer index @vector specifying units to be read. 
#    If @NULL, all units are read.}
#  \item{...}{Arguments passed to the low-level function for read units, 
#    e.g. @see "affxparser::readCelUnits" or @see "aroma.apd::readApdUnits".
#    If \code{units} is not already a CDF @list structure, these arguments
#    are also passed to \code{readUnits()} of @see "AffymetrixCdfFile".}
#  \item{force}{If @TRUE, cached values are ignored.}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a named @list structure.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getUnitIntensities", "AffymetrixCelSet", function(this, units=NULL, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="getUnitIntensities", class=class(this)[1], units=units, ...);
  id <- digest2(key);
  res <- this$.getUnitIntensitiesCache[[id]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "getUnitIntensitiesCache(): Returning cached data");
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the pathnames of all CEL files
  pathnames <- unlist(lapply(this, getPathname), use.names=FALSE);

  # Is 'units' a (pre-created) CDF structure?
  if (is.list(units)) {
    cdfUnits <- units;
  } else {
    # Always ask for CDF information from the CDF object!
    cdf <- getCdf(this);
    cdfUnits <- readUnits(cdf, units=units, ...);
  }

  res <- readCelUnits(pathnames, cdf=cdfUnits, readStdvs=FALSE, 
                              readPixels=FALSE, dropArrayDim=FALSE, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    this$.getUnitIntensitiesCache <- list();
    this$.getUnitIntensitiesCache[[id]] <- res;
    verbose && cat(verbose, "readUnits(): Updated cache");
  }

  res;
})



setMethodS3("readUnits", "AffymetrixCelSet", function(this, units=NULL, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "readCelUnits() of AffymetrixCelSet");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Generating hashcode key for cache");
  key <- list(method="readUnits", class=class(this)[1]);
  if (is.list(units)) {
    # BUG FIX: Cannot stratify only on unit names, e.g. then 'stratifyBy'
    # will not make a difference! 2008-03-11
#    key <- c(key, units=names(units), ...);
#   ...will also not work (stratifyBy="pm" and "mm" gives the same)!
#    key <- c(key, unitNames=names(units), cdfUnitSize=object.size(units), ...);
    # This will allow us to pass pre-digested object.
    unitsHashcode <- attr(units, "hashcode");
    if (is.null(unitsHashcode)) {
      unitsHashcode <- digest2(units);
      attr(units, "hashcode") <- unitsHashcode;
    }
    key <- c(key, unitsHashcode=unitsHashcode, ...);
  } else {
    key <- c(key, units=units, ...);
  }
  id <- digest2(key);
  verbose && exit(verbose);
  if (!force) {
    verbose && enter(verbose, "Trying to obtain cached data");
    res <- this$.readUnitsCache[[id]];
    verbose && exit(verbose);
    if (!is.null(res)) {
      verbose && cat(verbose, "readUnits(): Returning cached data");
      return(res);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the pathnames of all CEL files
  pathnames <- getPathnames(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read data from file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Calling readCelUnits() for ", 
                                              length(pathnames), " files");
  if (is.list(units)) {
    res <- readCelUnits(pathnames, cdf=units, dropArrayDim=FALSE, ...);
  } else {
    # Always ask for CDF information from the CDF object!
    verbose && enter(verbose, "Retrieving CDF unit information");
    cdf <- getCdf(this);
    cdfList <- getCellIndices(cdf, units=units, ..., verbose=less(verbose));
#    suppressWarnings({
#      cdfList <- readUnits(cdf, units=units, ..., verbose=less(verbose));
#    });
    verbose && str(verbose, cdfList[1]);
    verbose && exit(verbose);
    verbose && enter(verbose, "Retrieving CEL units across samples");
    res <- readCelUnits(pathnames, cdf=cdfList, dropArrayDim=FALSE, ...);
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    this$.readUnitsCache <- list();
    this$.readUnitsCache[[id]] <- res;
    verbose && cat(verbose, "readUnits(): Updated cache");
  }

  verbose && exit(verbose);

  res;
})


###########################################################################/**
# @RdocMethod getAverageFile
#
# @title "Calculates the mean and the standard deviation of the cell signal (intensity, standard deviation etc.) across the CEL set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{name}{The label of the calculated parameters.
#    If @NULL, a default name format \code{<prefix>-<mean>-<sd>} is used.}
#  \item{indices}{An @integer @vector specifying which cells to consider.
#    If \code{"remaining"}, only parameters for cells that have not been
#    are calculated.
#    If @NULL, all cells are used.}
#  \item{mean}{A @character of a @function specifying the function used
#    to calculate the average.}
#  \item{sd}{A @character of a @function specifying the function used
#    to calculate the standard deviation.}
#  \item{na.rm}{If @TRUE, @NAs are excluded before, otherwise not.}
#  \item{...}{Not used.}
#  \item{cellsPerChunk}{A @integer specifying the total number of cells 
#    (across arrays) read into memory per chunk.}
#  \item{moreCells}{A @double scalar indicating if more or less cells
#    should be used per chunk.}
#  \item{force}{If @TRUE, parameters for cells already calculated are
#    recalculated, otherwise not.}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns an @see "AffymetrixCelSet" of the same class as the CEL set
#   averaged.
# }
#
# \details{
#   The parameter estimates are stored as a CEL file of the same class as
#   the data files in the set.  The CEL file is named \code{<name>.cel}
#   and placed in the directory of the set.
#   Currently there is no specific data class for this file, but the average
#   cell signals are stored as "intensities", the standard deviation of the
#   cell signals as "stddevs", and the number of data points used for each
#   estimate is stored as "pixels".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAverageFile", "AffymetrixCelSet", function(this, name=NULL, prefix="average", indices="remaining", field=c("intensities", "stdvs"), mean=c("median", "mean"), sd=c("mad", "sd"), na.rm=TRUE, g=NULL, h=NULL, ..., cellsPerChunk=moreCells*10^7/length(this), moreCells=1, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'field':
  field <- match.arg(field);
  
  # Argument 'mean':
  if (is.character(mean)) {
    mean <- match.arg(mean);
    meanName <- mean;
    if (mean == "mean") {
      mean <- base::rowMeans;
    } else if (mean == "median") {
      mean <- rowMedians;
    }
  } else if (is.function(mean)) {
    meanName <- "customMean";
  } else {
    throw("Argument 'mean' must be either a character or a function: ", mode(mean));
  }

  # Argument 'sd':
  if (is.character(sd)) {
    sd <- match.arg(sd);
    sdName <- sd;
    if (sd == "sd") {
      sd <- rowSds;
    } else if (sd == "mad") {
      sd <- rowMads;
    }
  } else if (is.function(sd)) {
    sdName <- "customSd";
  } else {
    throw("Argument 'sd' must be either a character or a function: ", 
                                                           mode(sd));
  }

  # Argument 'name':
  if (is.null(name)) {
    key <- list(method="getAverageFile", class=class(this)[1], 
                arrays=sort(getNames(this)), mean=meanName, sd=sdName);
    # assign mean and sd to an empty environment so that digest() doesn't
    # pick up any "promised" objects from the original environment.
    # A bit ad hoc, but it works for now. /2007-01-03
    key <- base::lapply(key, FUN=function(x) {
      if (is.function(x))
        environment(x) <- emptyenv();
      x;
    })
    id <- digest2(key);
    name <- sprintf("%s-%s-%s-%s,%s", prefix, field, meanName, sdName, id);

    # Save for troubleshooting inconsistency(?) in digest()/Rv2.6.0?!?
    # /HB 2007-08-30
    idPath <- getCachePath(c("aroma.affymetrix", "idChecks"));
    idPathname <- file.path(idPath, id);
    saveObject(list(key=key, keyIds=lapply(key, digest2), id=id), idPathname);
  }

  # Argument 'indices':
  df <- as.list(this)[[1]];
  nbrOfCells <- getHeader(df)$total;
  if (force) {
    if (identical(indices, "remaining")) {
      indices <- NULL;
    }
  }

  if (is.null(indices)) {
    indices <- 1:nbrOfCells; 
  } else if (identical(indices, "remaining")) {
  } else {
    indices <- Arguments$getIndices(indices, range=c(1, nbrOfCells));
  }

  # Argument 'cellsPerChunk':
  cellsPerChunk <- Arguments$getInteger(cellsPerChunk, range=c(1,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving average cell signals across ", length(this), " arrays");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create CEL file to store the average array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create a private filename (with a dot prefix) to make sure it is not
  # identified as a regular CEL file when the directory is scanned for files.
  filename <- sprintf(".%s.CEL", name);
  if (is.null(this$.averageFiles))
    this$.averageFiles <- list();
  res <- this$.averageFiles[[filename]];

  # Has file been deleted since last time?
  if (!is.null(res) && !isFile(res)) {
    warning("Will recalculate average file, because it seems to have been deleted since last time: ", getPathname(res));
    res <- NULL;
  }

  if (is.null(res)) {
    verbose && enter(verbose, "Creating CEL file to store average signals");
    verbose && cat(verbose, "Pathname: ", file.path(getPath(this), filename));
    res <- createFrom(df, filename=filename, path=getPath(this), methods="create", clear=TRUE, verbose=less(verbose));
    verbose && exit(verbose);
    this$.averageFiles[[filename]] <- res;
  }

  verbose && print(verbose, res);

  pathname <- getPathname(res);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify which indices to use
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (identical(indices, "remaining")) {
    pixels <- readCel(pathname, readIntensities=FALSE, readStdvs=FALSE, 
                      readPixels=TRUE)$pixels;
    indices <- which(pixels == 0);
    rm(pixels); # Not needed anymore.
  }

  nbrOfIndices <- length(indices);

  # Nothing more to do?
  if (nbrOfIndices == 0)
    return(res);

  verbose && cat(verbose, "Number of cells to be updated: ", nbrOfIndices);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate the mean and standard deviation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Since we might want to do this robustly, but also because we want to
  # estimate the standard deviation, for each cell we need all data across 
  # arrays at once.  In order to this efficiently, we do this in chunks
  idxs <- 1:nbrOfIndices;
  head <- 1:cellsPerChunk;
  nbrOfChunks <- ceiling(nbrOfIndices / cellsPerChunk);
  verbose && cat(verbose, "Number cells per chunk: ", cellsPerChunk);

  # Get the pathnames of all CEL files to average
  pathnames <- lapply(this, getPathname);
  pathnames <- unlist(pathnames, use.names=FALSE);
  nbrOfArrays <- length(pathnames);

  if (!na.rm)
    n <- rep(nbrOfArrays, length=cellsPerChunk);
  count <- 1;
  while (length(idxs) > 0) {
    verbose && enter(verbose, "Fitting chunk #", count, " of ", nbrOfChunks);
    if (length(idxs) < cellsPerChunk) {
      head <- 1:length(idxs);
      if (!na.rm)
        n <- rep(nbrOfArrays, length=length(idxs));
    }

    # The indices to be used in this chunk
    ii <- idxs[head];
    verbose && cat(verbose, "Chunk size: ", length(ii));

    verbose && enter(verbose, "Reading data");
#    X <- readCelIntensities(pathnames, indices=indices[ii]);
    readIntensities <- (field == "intensities");
    readStdvs <- (field == "stdvs");
    # TODO: Ideally, affxparser::readCel() should support 
    # multiple filenames turning every data fields into a 
    # matrix. /HB 2007-01-07
    X <- matrix(as.double(NA), nrow=length(ii), ncol=nbrOfArrays);
    for (kk in seq(length=nbrOfArrays)) {
      X[,kk] <- readCel(filename = pathnames[kk],
                        indices = indices[ii],
                        readIntensities = readIntensities,
                        readHeader = FALSE,
                        readStdvs = readStdvs,
                        readPixels = FALSE,
                        readXY = FALSE,
                        readOutliers = FALSE,
                        readMasked = FALSE,
                        ...,
                        verbose = (verbose - 1))[[field]];
    }
    verbose && exit(verbose);

    if (!is.null(g)) {
      verbose && enter(verbose, "Transforming data using y = g(x)");
      X <- g(X);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Estimating averages and standard deviations");
    if (na.rm)
      n <- base::apply(X, MARGIN=1, FUN=function(x) { sum(!is.na(x)) });

    # Calculate the mean signal    
    mu <- mean(X, na.rm=na.rm);          # Special mean()!
    # Calculate the standard deviation of the signals
    sigma <- sd(X, mean=mu, na.rm=na.rm);   # Special sd()!

    verbose && exit(verbose);

    if (!is.null(h)) {
      verbose && enter(verbose, "Back-transforming estimates using x = h(y)");
      mu <- h(mu);
      sigma <- h(sigma);
      verbose && exit(verbose);
    }

    # Write estimates to result file
    verbose && enter(verbose, "Writing estimates");
    updateCel(pathname, indices=indices[ii], intensities=mu, stdvs=sigma, pixels=n);
    verbose && exit(verbose);

    # Not needed anymore
    mu <- sigma <- NULL;

    # Next chunk...
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collection
    gc();
    verbose && exit(verbose);
  } # while()

  verbose && exit(verbose);

  res;  
})


setMethodS3("getAverage", "AffymetrixCelSet", function(this, ...) {
  getAverageFile(this, ...);
})

setMethodS3("getAverageLog", "AffymetrixCelSet", function(this, ...) {
  getAverageFile(this, g=log2, h=function(x) 2^x, ...);
})

setMethodS3("getAverageAsinh", "AffymetrixCelSet", function(this, ...) {
  getAverageFile(this, g=asinh, h=sinh, ...);
})



setMethodS3("range", "AffymetrixCelSet", function(this, ...) {
  range(unlist(lapply(this, FUN=range, ...), use.names=FALSE));
})



setMethodS3("applyToUnitIntensities", "AffymetrixCelSet", function(this, units=NULL, FUN, stratifyBy="pm", verbose=FALSE, ...) {
  y <- getUnitIntensities(this, units=units, stratifyBy=stratifyBy, ...);

  y <- base::lapply(y, FUN=function(unit) {
    groups <- base::lapply(unit, FUN=function(group) {
      FUN(group[[1]], ...)
    })
    groups;
  })
  y;
}, private=TRUE)


setMethodS3("[", "AffymetrixCelSet", function(this, units=NULL, ..., drop=FALSE) {
  res <- readUnits(this, units=units, ...);
  if (drop && length(res) == 1)
    res <- res[[1]];
  res;
})

setMethodS3("[[", "AffymetrixCelSet", function(this, units=NULL, ...) {
  this[units=units, ..., drop=TRUE];
})


setMethodS3("getFullName", "AffymetrixCelSet", function(this, parent=1, ...) {
  NextMethod("getFullName", this, parent=parent, ...);
})


setMethodS3("getUnitGroupCellMap", "AffymetrixCelSet", function(this, ...) {
  ce <- getFile(this, 1);
  getUnitGroupCellMap(ce, ...);
})



############################################################################
# HISTORY:
# 2009-05-23
# o Now the chip type validation of fromFiles() for AffymetrixCelSet
#   is aware of tags in the chip type of the CEL files. This may happen
#   if custom CDFs are used are their full chip types are stored in the
#   CEL files, e.g. Hs_PromPR_v02,Harvard,ROIs,unique.
# 2008-12-18
# o BUG FIX: getUnitIntensities() of AffymetrixCelSet would drop the array
#   dimension if only one array was read.
# 2008-07-21
# o Now findByName() assert that the data set name is not empty.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: getData(), getAverageFile().
# 2008-07-14
# o Added protected update().
# o Added private updateSampleAnnotationSet().
# 2008-06-22
# o BUG FIX: Added missing getPlatform().
# 2008-05-31
# o BUG FIX: readUnits() would throw 'Error in readCelUnits(pathnames, cdf
#   = cdf, ...) : No CDF file for chip type found: GenomeWideSNP_6', if
#   the CDF was set to GenomeWideSNP_6,Full.CDF and no GenomeWideSNP_6.cdf 
#   file was found.  This was because readUnits() retrieve (x,y) information
#   units, which requires affxparser::readCelUnits() to locate the CDF to
#   infer the number of probe column in order to map (x,y) to cell indices.
#   Now readUnits() get the cell indices directly using getCellIndices().
#   Thanks Yue Hu for noticing this problem.
# 2008-05-17
# o Added getChipType().
# 2008-05-09
# o Now AffymetrixCelSet inherits from AromaMicroarrayDataSet.
# 2008-05-08
# o Now we are using foo(static, ...) instead of static$foo(...).
# o If paths=NULL in findByName(), it becomes the default argument value.
# o BUG FIX: fromFiles() was not declared static.
# o Made fromFiles() protected.
# 2008-03-11
# o BUG FIX: Calling readUnits(..., units=cdfUnits) twice, where 'cdfUnits' 
#   was first created with, say, stratifyBy="pm" and then with
#   stratifyBy="pmmm", for the identical set of units, would give the same
#   results both times.  This was because the hash key list object for the
#   file cache where identical.
# o BUG FIX: getUnitIntensities() would not pass arguments '...', e.g.
#   'stratifyBy', to readUnits() for AffymetrixCdfFile.  Thanks Tim 
#   Keighley, CSIRO, Sydney for reporting this.
# 2008-03-05
# o Updated message specifying that it looks for *sample* annotation files.
# 2008-02-28
# o Added getUnitGroupCellMap() to all AffymetrixCelSet classes.
# 2008-02-25
# o Now fromFiles() of AffymetrixCelSet uses the file cache, which speed up
#   the setup for extremely large data sets.  Note, currently the default
#   is 'force=TRUE' (for safety).
# 2008-01-30
# o Now as.character() for AffymetrixCelSet only reports time stamp if the
#   number of arrays in the data set is less that the threshold specified
#   in the 'aroma.affymetrix' settings.
# 2008-01-11
# o ROBUSTNESS: It was possible to set a non-compatible CDF when using
#   static fromFiles() of AffymetrixCelSet.
# o Added argument '.checkArgs' to setCdf() of AffymetrixCelSet, which now
#   validated the arguments once and not for every file.
# 2007-12-08
# o Added argument 'cdf' to static fromName() of AffymetrixCelSet.  When
#   using this argument, the 'chipType' argument is optional, and the 
#   returned CEL set will be using the specified CDF.
# o Now fromFiles() of AffymetrixCelSet takes argument 'cdf' which can be
#   used to override the default CDF.  This way the default CDF must not
#   have to be in the search path.
# 2007-09-25
# o Now getAverageFile() will detect if an average file has been deleted 
#   between calls and recalculate it.
# 2007-09-13
# o Added validate that the 'name' of an AffymetrixCelSet and 
#   AffymetrixCelFile is of at least length one.
# 2007-09-06
# o Now setCdf() throws an (informative) error message whenever one tries
#   to use an ASCII CDF file. This behavior can be changed by setting
#   rule 'rules$allowAsciiCdfs' in the package settings.
# 2007-08-01
# o Renamed static fromName() of AffymetrixCelSet to byName().
# 2007-04-06
# o BUG FIX: fromFiles() of AffymetrixCelSet would give error "Exception: 
#   Pathname not found: annotationData/samples" if that directory was 
#   missing.  Now it is instead created.
# 2007-04-03
# o Now fromFiles() verifies that the set CDF is compatible with the CEL
#   files, otherwise the CDF of the first CEL file is used.
# 2007-04-02
# o Now sample annotation files are searched for in annotationData/samples/.
# 2007-03-28
# o Added argument 'cache=TRUE' to getUnitIntensities() and readUnits().
# 2007-03-24
# o BUG FIX: clearCache() did not clear the .readUnitsCache field.
# 2007-03-16
# o BUG FIX: getAverageFile() of AffymetrixCelSet would average the wrong
#   set of cells if argument 'indices' was different from NULL.
# 2007-03-06
# o Now attributes are set from SAF files in fromFiles().
# 2007-02-22
# o Fixed the warning about "'tzone' attributes are inconsistent". See
#   code of as.character() for explanation.
# o Now fromFiles() accepts argument 'chipType' to override any chip type
#   specified in the CEL headers. This is useful in case different CEL files
#   refers to different chip types, which can be the case for mixed 
#   generations of CEL files.  Also added a scan of chip types.
# 2007-02-14
# o Added test for correct directory structure to fromFiles().  This will
#   enforce users to use the correct structure so that for instance the
#   name of the data set is correctly inferred.
# o Added argument 'onDuplicates' to fromFiles() so it is possible to 
#   exclude duplicated CEL files.
# 2007-02-06
# o Added static method findByName() and fromName().
# 2007-01-15
# o Added 'classes=class(this)' to all "digest" keys.
# 2007-01-07
# o BUG FIX: In KS's update of getAverageFile() to support averaging
#   over other fields than intensities, argument 'indices' was missing
#   in the readCel() call making the function fail when processed chunk
#   by chunk. /HB
# 2007-01-05
# o Removed getSampleNames().
# 2006-12-01
# o Now as.character() reports the range of CEL header timestamps.
# 2006-11-07
# o Now getAverageFile() uses rowMedians() of R.native if available, 
#   otherwise a local version utilizing apply(). Same for rowMads().
# 2006-10-24
# o Added getAverageLog() and getAverageAsinh().
# o Added transforms and anti-transforms g() and h() to getAverageFile().
# o Changed the defaults from mean to median, and sd to mad for 
#   getAverageFile().
# o Added Rdoc comments to getAverageFile().
# 2006-10-10
# o Renamed rma and gcrma to rmaSummary and gcrmaSummary, to avoid clash
#   with existing functions.
# o Added gcrma() wrapper function.
# o Added rma() wrapper function.
# o Fixed a bug in getData() - default for argument "fields" contained "xy",
#   which is not a valid field (x, y are separate).
# 2006-10-02
# o Added getData().  Now getIntensities() works again and is just a wrapper
#   to getData().
# 2006-09-18
# o Now references to all requested average files are cached so it can
#   return the same object instead of creating a new one each time.
# 2006-09-16
# o Added getSiblings() to easily get other data sets for the same
#   samples.
# 2006-09-14
# o Added a read-buffer cache to readUnits() and getUnitIntensities().
# 2006-08-27
# o Added getAverageFile().
# 2006-08-26
# o Now getName() of a CEL set is inferred from the pathname:
#     path/to/<name>/chip_files/<"chip type">/
# 2006-08-21
# o Now AffymetrixCelSet inherits from AffymetrixFileSet.
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
