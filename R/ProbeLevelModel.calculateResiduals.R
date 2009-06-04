setMethodS3("calculateResidualSet", "ProbeLevelModel", function(this, units=NULL, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  ces <- getChipEffectSet(this);
  if (inherits(ces, "CnChipEffectSet")) {
    if (ces$combineAlleles) {
      throw("calculateResidualSet() does not yet support chip effects for which allele A and allele B have been combined.");
    }
  }
  paf <- getProbeAffinityFile(this);
  nbrOfArrays <- nbrOfArrays(ces);

  # If residuals already calculated, and if force==FALSE, just return
  # a CelSet with the previous calculations

  verbose && enter(verbose, "Calculating PLM residuals");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get data and parameter objects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this);
  if (is.null(ds)) {
    throw("No data set specified for PLM: ", getFullName(this));
  }

  # Get the function how to calculate residuals.
  # Default is eps = y - yhat, but for instance RMA uses eps = y/yhat.
  calculateEps <- getCalculateResidualsFunction(this);

  cdf <- getCdf(ds);
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
  } else {
    nbrOfUnits <- length(units);
  }
  verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits);

  cdfData <- NULL;
  chipType <- getChipType(cdf);
  key <- list(method="calculateResiduals", class=class(this)[1], 
              chipType=chipType, params=getParameters(this),
              units=units);
  dirs <- c("aroma.affymetrix", chipType);
  if (!force) {
    cdfData <- loadCache(key, dirs=dirs);
    if (!is.null(cdfData)) {
      names(cdfData) <- gsub("cells2", "ceCells", names(cdfData), fixed=TRUE);
      verbose && cat(verbose, "Found indices cached on file");
    }
  }

  if (is.null(cdfData)) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identifying (unitGroupSizes, cells, ceCells)
    #
    # Note: This will take several minutes, but results will be save to 
    # file cache, so it is basically only done once.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    units0 <- units;
    if (is.null(units)) {
      units <- seq(length=nbrOfUnits(cdf));
    }
    nbrOfUnits <- length(units);

    # Memory optimization: Do things in chunks
    unitChunks <- splitInChunks(units, chunkSize=100e3);
    cdfData <- list(unitGroupSizes=NULL, cells=NULL, ceCells=NULL);
    for (kk in seq(along=unitChunks)) {
      verbose && enter(verbose, sprintf("Chunk #%d of %d", kk, length(unitChunks)));
      units <- unitChunks[[kk]];
      verbose && enter(verbose, "Retrieving CDF cell indices");
      cdfUnits <- getCellIndices(this, units=units, verbose=less(verbose));
      names(cdfUnits) <- NULL;  # Saves memory.
      gc <- gc();
      verbose && exit(verbose);
      verbose && enter(verbose, "Calculate group sizes");
      unitGroupSizes <- applyCdfGroups(cdfUnits, base::lapply, FUN=function(group) {
        length(.subset2(group, 1));
      });
      unitGroupSizes <- unlist(unitGroupSizes, use.names=FALSE);
      verbose && str(verbose, unitGroupSizes);

      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);

      cells <- unlist(cdfUnits, use.names=FALSE);

      rm(cdfUnits);
      # Garbage collect
      gc <- gc();

      verbose && enter(verbose, "Retrieving CDF cell indices for chip effects");
      cdfUnits <- getCellIndices(ces, units=units, verbose=less(verbose));
      ceCells <- unlist(cdfUnits, use.names=FALSE);
      verbose && exit(verbose);

      rm(cdfUnits); # Not needed anymore

      # Store
      cdfData$unitGroupSizes <- c(cdfData$unitGroupSizes, unitGroupSizes);
      cdfData$cells <- c(cdfData$cells, cells);
      cdfData$ceCells <- c(cdfData$ceCells, ceCells);

      verbose && str(verbose, cdfData);

      rm(unitGroupSizes, cells, ceCells);

      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);

      verbose && exit(verbose);
    } # for (kk in ...)

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    units <- units0;
    rm(units0);

    verbose && enter(verbose, "Saving to file cache");
    saveCache(cdfData, key=key, dirs=dirs);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  }

  verbose && cat(verbose, "CDF related data cached on file:");
  unitGroupSizes <- cdfData$unitGroupSizes;
  verbose && cat(verbose, "unitGroupSizes:");
  verbose && str(verbose, unitGroupSizes);
  cells <- cdfData$cells;
  verbose && cat(verbose, "cells:");
  verbose && str(verbose, cells);
  ceCells <- cdfData$ceCells;
  verbose && cat(verbose, "ceCells:");
  verbose && str(verbose, ceCells);
  rm(cdfData);
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Validate correctness
  if (!identical(length(unitGroupSizes), length(ceCells))) {
    throw("Internal error: 'unitGroupSizes' and 'ceCells' are of different lengths: ", length(unitGroupSizes), " != ", length(ceCells));
  }

  # Optimized reading order
  o <- .Internal(qsort(cells, TRUE));
  cells <- o$x;
  o <- o$ix;
  oinv <- .Internal(qsort(o, TRUE))$ix;
  
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate residuals for each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- getPath(this);
  phi <- NULL;
  for (kk in seq(ds)) {
    # Get probe-level signals
    df <- getFile(ds, kk);

    # Get chip effect estimates
    cef <- getFile(ces, kk);

    verbose && enter(verbose, sprintf("Array #%d ('%s')", kk, getName(df)));

    filename <- sprintf("%s,residuals.CEL", getFullName(df));
    pathname <- Arguments$getWritablePathname(filename, path=path);

    # Rename lower-case *.cel to *.CEL, if that is the case.  Old versions
    # of the package generated lower-case CEL files. /HB 2007-08-09
    pathname <- AffymetrixFile$renameToUpperCaseExt(pathname);

    verbose && cat(verbose, "Pathname: ", pathname);
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already calculated.");
      verbose && exit(verbose);
      next;
    }

    verbose && enter(verbose, "Retrieving probe intensity data");
    y <- getData(df, indices=cells, fields="intensities")$intensities[oinv];
    verbose && exit(verbose);

    if (is.null(phi)) {
      verbose && enter(verbose, "Retrieving probe-affinity estimates");
      phi <- getData(paf, indices=cells, fields="intensities")$intensities[oinv];
      verbose && exit(verbose);
    }

    # Assert that there is the correct number of signals 
    if (length(y) != length(phi)) {
      throw("Internal error: 'y' and 'phi' differ in lengths: ", 
                                           length(y), " != ", length(phi));
    }

    verbose && enter(verbose, "Retrieving chip-effect estimates");
    theta <- getData(cef, indices=ceCells, fields="intensities")$intensities;
    theta <- rep(theta, times=unitGroupSizes);
    verbose && exit(verbose);

    # Assert that there is the correct number of (expanded) chip effects
    if (length(theta) != length(phi)) {
      throw("Internal error: 'theta' and 'phi' differ in lengths: ", 
                                       length(theta), " != ", length(phi));
    }

    verbose && enter(verbose, "Calculating residuals");
    yhat <- phi * theta;

    # Assert that there is the correct number of "predicted" signals 
    if (length(yhat) != length(y)) {
      throw("Internal error: 'yhat' and 'y' differ in lengths: ", 
                                           length(yhat), " != ", length(y));
    }

    eps <- calculateEps(y, yhat);  # Model class specific.
    verbose && str(verbose, eps);

    # Assert that there is the correct number of residuals
    if (length(eps) != length(y)) {
      throw("Internal error: 'eps' and 'y' differ in lengths: ", 
                                           length(eps), " != ", length(y));
    }

    verbose && exit(verbose);
    rm(y, yhat, theta);

    # Memory optimization: One less copy of 'eps'.
    eps <- eps[o];
    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && enter(verbose, "Storing residuals");
    tryCatch({
      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating empty CEL file for results, if missing");
      createFrom(df, filename=pathname, path=NULL, 
                     methods="create", clear=TRUE, verbose=less(verbose));
      verbose && exit(verbose);

      verbose && enter(verbose, "Writing residuals");
      updateCel(pathname, indices=cells, intensities=eps);
      verbose && exit(verbose);
    }, interrupt = function(intr) {
      verbose && print(verbose, intr);
      file.remove(pathname);
    }, error = function(ex) {
      verbose && print(verbose, ex);
      file.remove(pathname);
    })
    verbose && exit(verbose);

#    verbose && enter(verbose, "Verifying");
#    eps2 <- getData(rf, indices=cells, fields="intensities")$intensities[oinv];
#    stopifnot(all.equal(eps, eps2, tolerance=.Machine$double.eps^0.25));
#    verbose && exit(verbose);
#    rm(eps2);

    rm(eps);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (kk ...)
  rm(cells, phi, unitGroupSizes);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Define residual set
  # Inherit the CDF from the input data set.
  cdf <- getCdf(ds);
  rs <- ResidualSet$fromFiles(path, cdf=cdf, ...);
  
  verbose && exit(verbose);

  invisible(rs);
})


setMethodS3("getCalculateResidualsFunction", "ProbeLevelModel", function(static, ...) {
  function(y, yhat) {
    y-yhat;
  }
}, static=TRUE, protected=TRUE)


setMethodS3("calculateResiduals", "ProbeLevelModel", function(this, ...) {
  calculateResidualSet(this, ...);
}, private=TRUE)



############################################################################
# HISTORY:
# 2007-12-08
# o Now calculateResidualSet() of ProbeLevelModel utilizes the new 'cdf'
#   argument of fromFiles() in AffymetrixCelSet.
# 2007-08-17
# o Now calculateResidualSet() of ProbeLevelModel only loads probe-affinity
#   estimates if needed, i.e. if residuals are already calculated this 
#   function will return faster now.
# o MEMORY OPTIMIZATION: calculateResidualSet() of ProbeLevelModel does
#   part of the work in chunks.
# 2007-08-09
# o calculateResidualSet() of ProbeLevelModel now creates CEL files with 
#   upper-case filename extension "*.CEL", not "*.cel".  The reason for this
#   is that some software don't recognize lower case filename extensions :(  
# 2007-04-12
# o BUG FIX: There was a if (TRUE) {} statement in calculateResidualSet() 
#   that was supposed to be if (!fource) {} in the release version.
# 2007-03-15
# o Renamed calculateResiduals() to calculateResidualSet().
# o BUG FIX: calculateResiduals() of ProbeLevelModel would give non-zero
#   residuals for cells not fitted by the PLM. 
# 2007-03-14
# o Optimized memory usage in calculateResiduals() further.
# 2007-03-06
# o Now calculateResiduals(), for which chip effects where allele A and B 
#   have been combined, gives an error, explaining that the feature is
#   still to be implemented.
# 2007-02-16
# o BUG FIX: Already calculated residuals would be recalculated and 
#   end up as an empty file.
# 2007-02-14 HB + KS
# o Now residuals can be calculated differently for different PLM classes.
#   This is done by overriding static getCalculateResidualsFunction().
# 2007-02-12 HB
# o Rewritten from KS:s code.
############################################################################
