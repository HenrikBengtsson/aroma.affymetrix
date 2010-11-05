###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod bgAdjustOptical
#
# @title "Applies optical background correction to a set of CEL files"
#
# \description{
#  @get "title".
#
#  Adapted from @see "gcrma::bg.adjust.optical" in the \pkg{gcrma} package.
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The location to save the adjusted data files.}
#   \item{minimum}{The minimum adjusted intensity.  Defaults to 1.}
#   \item{subsetToUpdate}{The indices of the probes to be updated.
#     If @NULL, all are updated.}
#   \item{typesToUpdate}{Types of probes to be updated.  For more details,
#     see argument \code{types} of \code{identifyCells()} for the
#     @see "AffymetrixCdfFile" class.}
#   \item{...}{Not used.}
#   \item{overwrite}{If @TRUE, already adjusted arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{.deprecated}{Internal argument.}
# }
#
# \value{
#  Returns the background adjusted @see "AffymetrixCelSet" object.
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#*/###########################################################################
setMethodS3("bgAdjustOptical", "AffymetrixCelSet", function(this, path=NULL, name="bgOptical", subsetToUpdate=NULL, typesToUpdate=NULL, minimum=1, overwrite=FALSE, skip=!overwrite, ..., verbose=FALSE, .deprecated=TRUE) {
  if (.deprecated) {
    throw("bgAdjustOptical() is deprecated.  Please use the OpticalBackgroundCorrection class");
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    # Path structure: /probeData/<dataSet,tags>/<chipType>/
    rootPath <- "probeData";
    chipType <- getChipType(cdf, fullname=FALSE);
    path <- file.path(rootPath, getFullName(this), chipType);
  }
  if (!is.null(path)) {
    # Verify this path (and create if missing)
    path <- Arguments$getWritablePath(path);
  }

  if (identical(getPath(this), path)) {
    throw("Cannot calibrate data file. Argument 'path' refers to the same path as the path of the data file to be calibrated: ", path);
  }
  mkdirs(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Identifying the probes to be updated");
  subsetToUpdate <- identifyCells(cdf, indices=subsetToUpdate,
                                                     types=typesToUpdate);
  verbose && exit(verbose);

  verbose && cat(verbose, "Adjusting for optical effect for ", length(subsetToUpdate), " probes");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # optical effect correction for each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  nbrOfArrays <- nbrOfArrays(this);
  verbose && enter(verbose, "Adjusting ", nbrOfArrays, " arrays");
  dataFiles <- list();
  for (kk in seq(this)) {
    verbose && enter(verbose, sprintf("Array #%d of %d", kk, nbrOfArrays));
    df <- getFile(this, kk);
    verbose && print(verbose, df);
    dataFiles[[kk]] <- bgAdjustOptical(df, path=path, subsetToUpdate=subsetToUpdate, typesToUpdate=NULL, minimum=minimum, verbose=less(verbose), .deprecated=.deprecated);
    verbose && exit(verbose);

    rm(df);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
  }
  verbose && exit(verbose);

  # CDF inheritance
  res <- newInstance(this, dataFiles);
  setCdf(res, getCdf(this));

  res;
}, private=TRUE)


###########################################################################/**
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
#   \item{path}{If an affinities vector is not specified,
#      gives the path to a file storing the affinities.}
# }
#
# \details{
#   This method is not constant in memory! /HB 2007-03-26
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#*/###########################################################################
setMethodS3("calculateParametersGsb", "AffymetrixCelSet", function(this, nbrOfPms=25000, affinities=NULL, path=NULL, ..., verbose=FALSE) {

  verbose <- Arguments$getVerbose(verbose);

  cdf <- getCdf(this);
  
  # get path to affinities  
  if (is.null(path)) {
    # try to find affinities file
    paths <- getPathnames(this)[1];
    paths <- getParent(paths);
    paths <- getParent(paths);
    paths <- paste(".",
                   paths,
                   "data/",
                   sep=";", collapse=";");

    pattern <- paste(getChipType(cdf, "-affinities.apa", sep=""));
    affinityFile <- findFiles(pattern=pattern, paths=paths, firstOnly=TRUE);
    if (is.null(affinityFile))
      throw("Could not locate probe affinities file: ", pattern);
  }

  verbose && enter(verbose, "Extracting PM indices");
  cells <- getCellIndices(cdf, useNames=FALSE, unlist=TRUE, verbose=less(verbose,2));
  pmCells <- cells[isPm(cdf, cache=FALSE)];
  rm(cells);
  pmCells <- sort(pmCells);
  verbose && exit(verbose);

  narray <- length(this);

  #  set.seed(1);
  #  was present in original gcrma code; left in here to allow for consistency
  #  check between old and new versions

  # get a sorted random subset of PM to use in parameter estimation
  pmCells.random <- sample(pmCells, size=nbrOfPms);
  pmCells.random <- sort(pmCells.random);
  rm(pmCells);   # Not needed anymore
  
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "Extracting ", nbrOfPms, " random PM intensities across CEL set");
  # make sure we don't just sample from a single array; avoids problems
  # if we happened to choose a low quality or otherwise aberrant array
  iarray <- sample(1:narray, size=nbrOfPms, replace=TRUE);

  # For each array, read the signals randomized for that array
  # Confirmed to give identical results. /HB 2007-03-26
  pathnames <- getPathnames(this);
  pm.random2 <- vector("double", nbrOfPms);
  for (aa in 1:narray) {
    verbose && enter(verbose, sprintf("Array #%d of %d", aa, narray));
    # Cells to be read for this array
    idxs <- which(iarray == aa);
    cells <- pmCells.random[idxs];
    pm.random2[idxs] <- readCel(pathnames[aa], indices=cells, 
                       readIntensities=TRUE, readStdvs=FALSE)$intensities;
    rm(idxs, cells);
    verbose && exit(verbose);
  }
  rm(iarray, pathnames);

##   pm.random22 <- pm.random2;
##   # TO DO: Don't read all probes from all arrays
##   pm.random <- readCelIntensities(getPathnames(this), indices=pmCells.random);
##   pm.random2 <- vector("double", nrow(pm.random));
##   for (i in 1:nbrOfPms) {
##     pm.random2[i] <- pm.random[i, iarray[i]];
##   }
##   verbose && exit(verbose);
##   # clean up
##   rm(pm.random, iarray);   # Not needed anymore
##   stopifnot(identical(pm.random2, pm.random22))

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting probe affinities and fitting linear model")

  if (is.null(affinities)) {
    aff <- readApd(affinityFile, indices=pmCells.random)$affinities;
  } else {
    aff <- affinities[pmCells.random];
  }
  rm(pmCells.random);  # Not needed anymore

  # Work on the log2 scale
  pm.random2 <- log2(pm.random2);  # Minimize memory usage.
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "Fitting the GCRMA background linear model");
  verbose && str(verbose, pm.random2);
  verbose && str(verbose, aff);
  fit1 <- lm(pm.random2 ~ aff);
  verbose && print(verbose, fit1);
  verbose && exit(verbose);

  verbose && exit(verbose);
  
  fit1$coef;
}, private=TRUE)



###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod bgAdjustGcrma
#
# @title "Applies probe sequence based background correction to a set of
# CEL files"
#
# \description{
#  @get "title".
#
#  Adapted from @see "gcrma::bg.adjust.gcrma" in the \pkg{gcrma} package.
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The path where to save the adjusted data files.}
#   \item{name}{Name of the set containing the background corrected files.}
#   \item{type}{The type of background correction.  Currently accepted types
#       are "fullmodel" (the default, uses MMs) and "affinities" (uses
#       probe sequence only).}
#   \item{indicesNegativeControl}{Locations of any negative control
#       probes (e.g., the anti-genomic controls on the human exon array).
#       If @NULL and type=="affinities", MMs are used as the negative
#       controls.}
#   \item{opticalAdjust}{If @TRUE, apply correction for optical effect,
#       as in @see "gcrma::bg.adjust.optical".}
#   \item{gsbAdjust}{Should we adjust for specific binding (defaults to
#        @TRUE)?}
#   \item{k}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{rho}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{stretch}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{fast}{If @TRUE, an ad hoc transformation of the PM is performed
#       (\code{gcrma::gcrma.bg.transformation.fast}).}
#   \item{overwrite}{If @TRUE, already adjusted arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{.deprecated}{Internal argument.}
# }
#
# \value{
#  Returns the background adjusted @see "AffymetrixCelFile" object.
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#
# \seealso{
#  @see "gcrma::bg.adjust.gcrma"
#  @seeclass
# }
#*/###########################################################################
setMethodS3("bgAdjustGcrma", "AffymetrixCelSet", function(this, path=NULL, name="bgGcrma", probePath=NULL, affinities=NULL, type="fullmodel",  indicesNegativeControl=NULL, opticalAdjust=TRUE, gsbAdjust=TRUE, k=6 * fast + 0.5 * (1 - fast), rho=0.7, stretch=1.15*fast + (1-fast), fast=TRUE, overwrite=FALSE, skip=!overwrite, ..., verbose=FALSE, .deprecated=TRUE) {
  if (.deprecated) {
    throw("bgAdjustGcrma() is deprecated.  Please use the GcRmaBackgroundCorrection class");
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    rootPath <- "probeData";
    chipType <- getChipType(cdf);
    path <- file.path(rootPath, getFullName(this), chipType);
  }
  if (!is.null(path)) {
    # Verify this path (and create if missing)
    path <- Arguments$getWritablePath(path);
  }

  if (identical(getPath(this), path)) {
    throw("Cannot calibrate data file. Argument 'path' refers to the same path as the path of the data file to be calibrated: ", path);
  }
  mkdirs(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate probe affinities, if not already existing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(affinities)) {
    verbose && enter(verbose, "Computing probe affinities (independent of data)");

    # Alternative 1: Using ACS annotation file
    affinities <- NULL;
    tryCatch({
      affinities <- computeAffinitiesByACS(cdf, ..., verbose=less(verbose));
    }, error = function(ex) {});

    if (is.null(affinities)) {
      # Alternative 2: Using old probe-tab files (deprecated)
      filename <- paste(getChipType(cdf), "-affinities.apa", sep="");
      pathname <- filePath(path, filename, expandLinks="any");
  
      if (isFile(pathname)) {
        verbose && enter(verbose, "Reading saved affinities: ", pathname);
        affinities <- readApd(pathname)$affinities;
        verbose && exit(verbose);
      } else {
        affinities <- computeAffinities(cdf, paths=probePath, ..., verbose=less(verbose));
        verbose && cat(verbose, "Saving affinities: ", pathname);
        writeApd(pathname, data=affinities, name="affinities");
      }
    }

    verbose && printf(verbose, "RAM: %.2fMB\n", 
                                         object.size(affinities)/1024^2);
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # optical background correction
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (opticalAdjust) {
    OBG <- OpticalBackgroundCorrection(this);
    dsOBG <- process(OBG, ..., verbose=verbose);
    rm(OBG);
    this <- dsOBG;
  }

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # estimate specific binding (GSB, in gcrma terminology)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  if (gsbAdjust) {
    verbose && enter(verbose, "Estimating specific binding parameters (data dependent)");
    parametersGsb <- calculateParametersGsb(this, affinities=affinities, path=path, ..., verbose=verbose);
    verbose && cat(verbose, "parametersGsb:");
    verbose && print(verbose, parametersGsb);
    verbose && exit(verbose);
  }
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # NSB correction for each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  nbrOfArrays <- nbrOfArrays(this);
  verbose && enter(verbose, "Adjusting ", nbrOfArrays, " arrays");
  dataFiles <- list();
  for (kk in seq(this)) {
    verbose && enter(verbose, sprintf("Array #%d of %d", kk, nbrOfArrays));
    df <- getFile(this, kk);
    verbose && print(verbose, df);
    dataFiles[[kk]] <- bgAdjustGcrma(df, path=path, type=type, indicesNegativeControl=indicesNegativeControl, affinities=affinities, gsbAdjust=gsbAdjust, parametersGsb=parametersGsb, k=k, rho=rho, stretch=stretch, fast=fast, overwrite=overwrite, skip=skip, ..., verbose=less(verbose), .deprecated=.deprecated);

    rm(df);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  res <- newInstance(this, dataFiles);
  setCdf(res, getCdf(this));

  res;
}, private=TRUE)



###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod bgAdjustRma
#
# @title "Applies RMA background correction to a set of
# CEL files"
#
# \description{
#  @get "title".
#
#  Adapted from @see "affy::bg.adjust" in the \pkg{affy} package.
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The path where to save the adjusted data files.}
#   \item{name}{Name of the set containing the background corrected files.}
#   \item{overwrite}{If @TRUE, already adjusted arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{...}{Not used.}
#   \item{.deprecated}{Internal argument.}
# }
#
# \value{
#  Returns the background adjusted @see "AffymetrixCelFile" object.
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#
# \seealso{
#  @see "affy::bg.adjust"
#  @seeclass
# }
#*/###########################################################################
setMethodS3("bgAdjustRma", "AffymetrixCelSet", function(this, path=NULL, tags="RBC", pmonly=TRUE, addJitter=FALSE, jitterSd=0.2, overwrite=FALSE, skip=!overwrite, ..., verbose=FALSE, .deprecated=TRUE) {
  if (.deprecated) {
    throw("bgAdjustRma() is deprecated.  Please use the RmaBackgroundCorrection class");
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    # Path structure: /bgRma/<data set name>/chip_data/<chip type>/
    rootPath <- "probeData";
    fullname <- paste(c(getFullName(this), tags), collapse=",");
    chipType <- getChipType(cdf);
    path <- file.path(rootPath, fullname, chipType);
  }
  if (!is.null(path)) {
    # Verify this path (and create if missing)
    path <- Arguments$getWritablePath(path);
  }

  if (identical(getPath(this), path)) {
    throw("Cannot calibrate data file. Argument 'path' refers to the same path as the path of the data file to be calibrated: ", path);
  }
  mkdirs(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # apply normal+exponential model to each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  nbrOfArrays <- nbrOfArrays(this);
  verbose && enter(verbose, "Adjusting ", nbrOfArrays, " arrays");
  dataFiles <- list();
  for (kk in seq(this)) {
    verbose && enter(verbose, sprintf("Array #%d of %d", kk, nbrOfArrays));
    df <- getFile(this, kk);
    verbose && print(verbose, df);
    dataFiles[[kk]] <- bgAdjustRma(df, path=path, pmonly=pmonly, addJitter=addJitter, jitterSd=jitterSd, overwrite=overwrite, skip=skip, ..., verbose=less(verbose), .deprecated=.deprecated);

    rm(df);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  res <- newInstance(this, dataFiles);
  setCdf(res, getCdf(this));

  res;
}, private=TRUE)


############################################################################
# HISTORY:
# 2010-09-29 [HB]
# o ROBUSTNESS: Now bgAdjustGcrma(..., affinities=NULL) is deprecated and
#   throws an exception.
# 2009-08-31 [HB]
# o CLEAN UP: Updated how 'path' is set internally if not specified.
# 2009-04-06 [HB]
# o BUG FIX: The output path of bgAdjustRma() of AffymetrixCelSet would
#   include the full chip type.  It would also give an error if no tags
#   where specified.
# 2009-03-29 [MR]
# o Made slight modifications for bgAdjustGcRma() to work with the 
#   newer Gene 1.0 ST arrays.
# 2007-09-06
# o Made calculateParametersGsb() more memory efficient, because it's using
#   the new unlist feature in getCellIndices() of AffymetrixCdfFile.
# 2007-06-30
# o Added .deprecated=TRUE to all methods.
# 2007-03-26
# o Speed up: Using isPm(cdf) instead of readCdfCellIndices() etc.
# o Memory optimization: Found an unlist() without use.names=FALSE.  Saves
#   about 95% of the object, or 250MB of RAM.
# o TO DO: Store GCRMA affinities under annotationData/ and not use the
#   data set paths. /HB
# o Tried to make calculateParametersGsb() constant in memory. Before it 
#   read all data for all arrays and then subsampled! /HB
# o BUG FIX: Argument 'minimum' was not used but forced to equal 1.
# 2007-03-22
# o rename gsbParameters to parametersGsb to avoid clash of arguments
#   in bgAdjustGcrma.AffymetrixCelFile().  Not sure why gsbAdjust and
#   gsbParameters are being matched, but there you go.
# 2006-10-10
# o add RMA background correction (normal+exponential)
# 2006-10-06
# o make sure cdf association is inherited
# 2006-10-04
# o Tested, debugged, docs added.
# 2006-09-28
# o Created (based on AffymetrixCelSet.NORM.R).
############################################################################
