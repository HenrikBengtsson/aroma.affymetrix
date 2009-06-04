###########################################################################/**
# @RdocClass AffymetrixCnChpSet
#
# @title "The AffymetrixCnChpSet class"
#
# \description{
#  @classhierarchy
#
#  A AffymetrixCnChpSet object represents a set of AffymetrixCnChpFile:s
#  with \emph{identical} chip types.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{files}{A @list of @see "AffymetrixCnChpFile":s.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \seealso{
#   @see "AffymetrixCnChpFile".
# }
#
# @author
#*/###########################################################################
setConstructorS3("AffymetrixCnChpSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    reqFileClass <- "AffymetrixCnChpFile";
    lapply(files, FUN=function(df) {
      if (!inherits(df, reqFileClass))
        throw("Argument 'files' contains a non-", reqFileClass, 
                                                  " object: ", class(df)[1]);
    })
  } else if (inherits(files, "AffymetrixCnChpSet")) {
    return(as.AffymetrixCnChpSet(files));
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }


  extend(AffymetrixFileSet(files=files, ...), "AffymetrixCnChpSet");
})



###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the set"
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
setMethodS3("as.character", "AffymetrixCnChpSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  s <- c(s, sprintf("Tags: %s", tags));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  n <- nbrOfArrays(this);
  s <- c(s, sprintf("Number of arrays: %d", n));
  names <- getNames(this);
  if (n >= 5)
    names <- c(names[1:2], "...", names[n]);
  names <- paste(names, collapse=", ");
  s <- c(s, sprintf("Names: %s", names));
  s <- c(s, sprintf("Total file size: %.2fMB", getFileSize(this)/1024^2));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)



setMethodS3("findByName", "AffymetrixCnChpSet", function(static, name, tags=NULL, chipType=NULL, paths=c("chpData"), ...) {
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
    paths <- eval(formals(findByName.AffymetrixCnChpSet)[["paths"]]);
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


setMethodS3("byName", "AffymetrixCnChpSet", function(static, name, tags=NULL, chipType=NULL, cdf=NULL, paths=NULL, ...) {
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


setMethodS3("fromFiles", "AffymetrixCnChpSet", function(static, path="rawData/", pattern="[.](cnchp|CNCHP)$", cdf=NULL, checkChipType=is.null(cdf), ..., fileClass="AffymetrixCnChpFile", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  
  verbose && enter(verbose, "Defining ", class(static)[1], " from files");

  set <- fromFiles.AffymetrixFileSet(static, path=path, pattern=pattern, ..., fileClass=fileClass, verbose=less(verbose));

  verbose && enter(verbose, "Retrieved files: ", nbrOfFiles(set));

  if (nbrOfFiles(set) > 0) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Scan all CHP files for possible chip types
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Chip type according to the directory structure
    path <- getPath(set);
    chipType <- basename(path);
    verbose && cat(verbose, 
                   "The chip type according to the path is: ", chipType);
  
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Use the same CDF object for all CEL files.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if (is.null(cdf)) {
      verbose && enter(verbose, "Retrieving the CDF for chip type '", chipType, "' inferred from path");
      cdf <- AffymetrixCdfFile$byChipType(chipType);
      verbose && exit(verbose);
  
      verbose && enter(verbose, "Check compatibility with 1st CEL file");
      verbose && cat(verbose, "Chip type: ", chipType);
      cf <- getFile(set, 1);
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

  verbose && exit(verbose);

  set;
})




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
setMethodS3("nbrOfArrays", "AffymetrixCnChpSet", function(this, ...) {
  nbrOfFiles(this, ...);
})



###########################################################################/**
# @RdocMethod as.AffymetrixCnChpSet
# @alias as.AffymetrixCnChpSet.list
# @alias as.AffymetrixCnChpSet.default
#
# @title "Coerce an object to an AffymetrixCnChpSet object"
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
#   Returns an @see "AffymetrixCnChpSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.AffymetrixCnChpSet", "AffymetrixCnChpSet", function(object, ...) {
  object;
})

setMethodS3("as.AffymetrixCnChpSet", "list", function(object, ...) {
  AffymetrixCnChpSet(object, ...);
})

setMethodS3("as.AffymetrixCnChpSet", "default", function(object, ...) {
  throw("Cannot coerce object to an AffymetrixCnChpSet object: ", mode(object));
})



setMethodS3("getFullName", "AffymetrixCnChpSet", function(this, parent=1, ...) {
  NextMethod("getFullName", this, parent=parent, ...);
})


setMethodS3("extractLogRatios", "AffymetrixCnChpSet", function(this, units=NULL, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'units':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract log ratios
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  readMap <- NULL;
  data <- NULL;
  nbrOfArrays <- nbrOfArrays(this);
  gcCount <- 0;
  for (kk in seq(length=nbrOfArrays)) {
    df <- getFile(this, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, getName(df), nbrOfArrays));

    if (!is.null(readMap))
      setUnitReadMap(df, readMap=readMap);

    dataKK <- extractLogRatios(df, units=units, ..., verbose=less(verbose, 5));

    if (is.null(readMap))
      readMap <- getUnitReadMap(df);

    verbose && str(verbose, dataKK);
    if (is.null(data)) {
      dim <- c(length(dataKK), nbrOfArrays);
      dimnames <- list(NULL, getNames(this));
      naValue <- as.double(NA);
      data <- array(naValue, dim=dim, dimnames=dimnames);
    }
    data[,kk] <- dataKK;
    rm(dataKK);

    # Garbage collect?
    gcCount <- gcCount + 1;
    if (gcCount %% 10 == 0) {
      gc <- gc();
      verbose && print(verbose, gc);
    }

    verbose && exit(verbose);
  } # for (kk ...)

  # Drop singleton dimensions
  if (drop) {
    data <- drop(data);
  }

  verbose && cat(verbose, "Log ratios:");
  verbose && str(verbose, data);

  data;
})


setMethodS3("getCdf", "AffymetrixCnChpSet", function(this, ...) {
  getCdf(this$files[[1]], ...);
})
 
setMethodS3("setCdf", "AffymetrixCnChpSet", function(this, cdf, verbose=FALSE, ..., .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    # Argument 'cdf':
    if (!inherits(cdf, "AffymetrixCdfFile")) {
      throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
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
 

############################################################################
# HISTORY:
# 2008-08-22
# o Created.
############################################################################
