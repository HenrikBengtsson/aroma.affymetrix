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
# @author
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
  probeModel <- match.arg(probeModel);

  extend(ParameterCelSet(...), "ChipEffectSet",
    "cached:.firstCells" = NULL,
    probeModel = probeModel
  )
})


setMethodS3("clearCache", "ChipEffectSet", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".firstCells")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)


setMethodS3("as.character", "ChipEffectSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod(generic="as.character", object=this, ...);
  params <- paste(getParametersAsString(this), collapse=", ");
  s <- c(s, sprintf("Parameters: (%s)", params));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getParameters", "ChipEffectSet", function(this, ...) {
  ce <- getFile(this, 1);
  getParameters(ce, ...);
})

setMethodS3("getParametersAsString", "ChipEffectSet", function(this, ...) {
  params <- getParameters(this);
  params <- trim(capture.output(str(params)))[-1];
  params <- gsub("^[$][ ]*", "", params);
  params <- gsub(" [ ]*", " ", params);
  params <- gsub("[ ]*:", ":", params);
  params;
}, private=TRUE)


setMethodS3("getChipEffectFileClass", "ChipEffectSet", function(static, ...) {
  ChipEffectFile;
}, static=TRUE, private=TRUE)


setMethodS3("findByName", "ChipEffectSet", function(static, ..., paths="plmData/") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'paths':
  if (is.null(paths)) {
    paths <- eval(formals(findByName.ChipEffectSet)[["paths"]]);
  }


  # Unfortunately method dispatching does not work here.
  path <- findByName.AffymetrixCelSet(static, ..., paths=paths);
  
  path;
}, static=TRUE)

setMethodS3("byPath", "ChipEffectSet", function(static, path="plmData/", pattern=",chipEffects[.](c|C)(e|E)(l|L)$", checkChipType=FALSE, ..., fileClass=NULL) {
  # Argument 'fileClass':
  if (is.null(fileClass))
    fileClass <- gsub("Set$", "File", class(static)[1]);

  # Unfortunately, method dispatching does not work here.
  byPath.AffymetrixCelSet(static, path=path, pattern=pattern, ..., fileClass=fileClass, checkChipType=checkChipType);
}, protected=TRUE, static=TRUE)


setMethodS3("fromDataSet", "ChipEffectSet", function(static, dataSet, path, name=getName(dataSet), cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Get the ChipEffectFile class specific for this set
  clazz <- getChipEffectFileClass(static);

  verbose && enter(verbose, "Retrieving chip-effects from data set");
  ces <- vector("list", length(dataSet));
  verbose && cat(verbose, "Data set: ", name);
  for (kk in seq(dataSet)) {
    df <- getFile(dataSet, kk);
    verbose && enter(verbose, 
                           sprintf("Retrieving chip-effect #%d of %d (%s)",
                                               kk, length(ces), getName(df)));
    ce <- clazz$fromDataFile(df, path=path, name=name, cdf=cdf, ..., 
                                                       verbose=less(verbose));
    if (is.null(cdf)) {
      verbose && enter(verbose, "Retrieving the CDF for the chip-effect file");
      cdf <- getCdf(ce);
      verbose && exit(verbose);
    }
    ces[[kk]] <- ce;
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Create an ChipEffectSet
  newInstance(static, ces);
})


setMethodS3("getCellIndices", "ChipEffectSet", function(this, ...) {
  # Use the first chip-effect file to get the CDF structure.
  # Note: Ideally we want to define a special CDF class doing this
  # instead of letting the data file do this. /HB 2006-12-18
  ce <- getFile(this, 1);
  getCellIndices(ce, ...);
})


setMethodS3("readUnits", "ChipEffectSet", function(this, units=NULL, cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading chip effects unit by unit for ", nbrOfArrays(this), " arrays");

  if (is.null(cdf)) {
    verbose && enter(verbose, "Getting cell indices from CDF");
    cdf <- getCellIndices(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  verbose && enter(verbose, "Calling readUnits() in superclass");
  res <- NextMethod("readUnits", this, units=cdf, ..., verbose=less(verbose));
  verbose && exit(verbose);

  # Get first chip-effect file and use that to decode the read structure
  # This takes some time for a large number of units /HB 2006-10-04
  ce <- getFile(this, 1);
  res <- decode(ce, res, verbose=less(verbose));

  verbose && exit(verbose);

  res;
})


setMethodS3("updateUnits", "ChipEffectSet", function(this, units=NULL, cdf=NULL, data, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Updating chip-effect files");

  # Get the CDF structure for all chip-effect files
  if (is.null(cdf)) {
    cdf <- getCellIndices(this, units=units, verbose=less(verbose, 1));
  }

  # Update each file one by one
  arrays <- seq(this);
  nbrOfArrays <- length(arrays);
  verbose && cat(verbose, "Number of files: ", nbrOfArrays);

  names <- getNames(this);

  verbose && enter(verbose, "Making sure the files are updated in lexicographic order");
  # Reorder such that the file with the "last" name is saved last
  fullnames <- getFullNames(this);
  o <- order(fullnames, decreasing=FALSE);
  arrays <- arrays[o];
  verbose && str(verbose, arrays);
  verbose && cat(verbose, "Last array: ", fullnames[arrays[nbrOfArrays]]);
  rm(fullnames, o);
  verbose && exit(verbose);

  verbose <- less(verbose);
  for (ii in arrays) {
    verbose && enter(verbose, sprintf("Array #%d of %d ('%s')",
                                         ii, nbrOfArrays, names[ii]));
    ce <- getFile(this, ii);

    verbose <- less(verbose, 50);
    verbose && enter(verbose, "Extracting estimates");  # 3-4s
    dataOne <- base::lapply(data, FUN=base::lapply, function(group) {
      # theta = group$theta[ii] = ...
      # stdvs = group$sdTheta[ii] = ...
      list(
        theta=.subset(.subset2(group, "theta"), ii),
        sdTheta=.subset(.subset2(group, "sdTheta"), ii),
        thetaOutliers=.subset(.subset2(group, "thetaOutliers"), ii)
      );
    });
    verbose && str(verbose, dataOne[1]);
    verbose && exit(verbose);

    # Assert correct lengths  /HB 2007-04-12
    if (length(cdf) != length(dataOne)) {
      throw("Internal error: Lengths of 'cdf' and 'dataOne' differ: ",
                                    length(cdf), " != ", length(dataOne));
    }

    verbose && enter(verbose, "Updating file");  # 6-7s ~98% in encode()
#    verbose && printf(verbose, "class(ce)[1]: %s\n", class(ce)[1]);
#    updateUnits(ce, cdf=cdf, data=dataOne, verbose=less(verbose, 50));
    updateUnits(ce, cdf=cdf, data=dataOne, verbose=verbose);
    rm(dataOne, ce);
    verbose && exit(verbose);

    verbose <- more(verbose, 50);

    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (ii ...)
  verbose <- more(verbose);
  verbose && exit(verbose);
}, protected=TRUE)



setMethodS3("getAverageFile", "ChipEffectSet", function(this, indices="remaining", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'indices':
  if (identical(indices, "remaining")) {
  } else if (is.null(indices)) {
#    # Update only cells which stores values
#    indices <- getCellIndices(this, verbose=verbose);
#    indices <- unlist(indices, use.names=FALSE);
  }

  NextMethod(generic="getAverageFile", object=this, indices=indices, ...);
})


setMethodS3("findUnitsTodo", "ChipEffectSet", function(this, ...) {
  # Look into the chip-effect file that comes last in a lexicographic
  # order, becuase that is updated last.
  names <- getFullNames(this);
  idx <- order(names, decreasing=TRUE)[1];
  df <- getFile(this, idx);
  findUnitsTodo(df, ...);
})



setMethodS3("extractMatrix", "ChipEffectSet", function(this, ..., field=c("theta", "sdTheta", "RLE", "NUSE"), drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'field':
  field <- match.arg(field);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  if (toupper(field) == "RLE") {
    verbose && enter(verbose, "Extracting chip effects and calculating RLE scores");

    # RLE - Relative Log2 Expression values
    # Get chip effect estimates (on the log scale)
    data <- extractMatrix(this, ..., field="theta", verbose=less(verbose, 1));
    data <- log2(data);  # ...stored on the intensity scale

    # Robust average (on the log scale, but stored on the intensity scale!)
    avg <- getAverageLog(this, field="intensities", mean="median", 
                                                      verbose=less(verbose,1));
    dataR <- extractMatrix(avg, field="theta", ..., verbose=less(verbose, 1));
    dataR <- log2(as.vector(dataR));

    # Log ratios of chip effects
    data <- data - dataR;
    rm(dataR, avg);

    verbose && exit(verbose);
  } else if (toupper(field) == "NUSE") {
    verbose && enter(verbose, "Extracting standard errors for chip effects and calculates NUSE scores");
    # NUSE - Normalized Unscaled Standard Errors
    # Get standard errors (on the log scale)
    data <- extractMatrix(this, ..., field="sdTheta", 
                                                   verbose=less(verbose, 1));
    data <- log2(data);  # ...stored on the intensity scale

    # Robust average (on the log scale, but stored on the intensity scale!)
    avg <- getAverageLog(this, field="stdvs", mean="median", 
                                                    verbose=less(verbose, 1));
    dataR <- extractMatrix(avg, field="theta", ..., verbose=less(verbose, 1));
    dataR <- log2(as.vector(dataR));

    # Log ratios of standard errors
    data <- data / dataR;
    rm(dataR, avg);

    verbose && exit(verbose);
  } else {
    data <- NextMethod("extractMatrix", this, ..., field=field);
  }

  # Drop singleton dimensions?
  if (drop) {
    data <- drop(data);
  }

  data;
})


############################################################################
# HISTORY:
# 2010-05-08
# o Now all findUnitsTodo() for data sets checks the data file that comes
#   last in a lexicographic ordering.  This is now consistent with how
#   the summarization methods updates the files.  Before it was use to be
#   the one that is last in the data set.
# 2010-02-20
# o ROBUSTNESS: Now updateUnits() of ChipEffectSet updates the files in
#   lexicographic order.  Before there was a risk that this was not done
#   if fullname translators are changing the lexicographic ordering.
# o MEMORY OPTIMIZATION: Now updateUnits() of ChipEffectSet cleans out the
#   temporary data object extracted for each chip-effect file written.
#   It also calls the garbage collector after each file written.
# 2008-07-09
# o Added argument drop=FALSE to extractMatrix().
# 2008-05-08
# o If paths=NULL in findByName(), it becomes the default argument value.
# o Made fromFiles() protected.
# 2008-02-25
# o Now extractMatrix() accepts special fields "RLE" and "NUSE".
# 2008-02-22
# o Now ChipEffectSet inherits from ParameterCelSet instead of as before
#   directly from AffymetrixCelSet.
# o Added extractMatrix() wrapper.
# 2007-12-10
# o Now fromDataSet() of ChipEffectSet accepts argument 'cdf'.
# 2007-07-01
# o BUG FIX: getOutputDataSet() of Transform would give "Error in 
#   fromFiles.AffymetrixCelSet(static, path = path, pattern = pattern,:
#   formal argument "checkChipType" matched by multiple actual arguments".
#   This was due to the recent adding of 'checkChipType=FALSE'.  Fixed
#   by adding 'checkChipType=FALSE' to fromFiles() of ChipEffectSet.
#   Thanks Jeremy Silver at WEHI for report and troubleshooting this.
# 2007-04-03
# o BUG FIX: Static fromFiles() did not call ditto in the super class but
#   instead in the grand-parent super class.
# 2007-02-19
# o Added findByName() and fromName().
# 2006-11-22
# o Updated fromFiles() so it automagically finds the file class.
# 2006-10-??
# o Updated fromFiles() to search for filename with tags 'chipEffects'.
#   Before files with suffix "-chipEffects' was searched for.
# 2006-09-10
# o Added findUnitsTodo().
# o Starting to make use of specially design CDFs and CEL files for storing
#   chip effects.  This make getFirstCellIndices() obsolete.
# 2006-08-28
# o Added getAverageFile() so that only cells that store actual chip-effect
#   estimates are averaged.
# 2006-08-26
# o Created.
############################################################################
