###########################################################################/**
# @set "class=AffymetrixCelFile"
# @RdocMethod bgAdjustOptical
#
# @title "Applies optical background correction to a CEL file"
#
# \description{
#  @get "title".
#
# Adapted from @see "gcrma::bg.adjust.optical" in the \pkg{gcrma} package.
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
#*/###########################################################################
setMethodS3("bgAdjustOptical", "AffymetrixCelFile", function(this, path=file.path("bgOptical", getChipType(this)), minimum=1, subsetToUpdate=NULL, typesToUpdate=NULL, overwrite=FALSE, skip=!overwrite, verbose=FALSE, ..., .deprecated=TRUE) {
  if (.deprecated) {
    throw("bgAdjustOptical() is deprecated.  Please use the OpticalBackgroundCorrection class");
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePathname(path);
  if (identical(getPath(this), path)) {
    throw("Cannot background correct data file. Argument 'path' refers to the same path as the path of the data file to be corrected: ", path);
  }


  cdf <- getCdf(this);
  nbrOfCells <- nbrOfCells(cdf);

  # Argument 'subsetToUpdate':
  getFraction <- (length(subsetToUpdate) == 1) &&
                               (subsetToUpdate >= 0) && (subsetToUpdate < 1);
  if (!getFraction) {
    subsetToUpdate <- Arguments$getIndices(subsetToUpdate,
                                           range=c(1, nbrOfCells));
  }
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- basename(getPathname(this));
  filename <- gsub("[.]cel$", ".CEL", filename);  # Only output upper case!
  pathname <- Arguments$getWritablePathname(filename, path=path,
                                         mustNotExist=(!overwrite && !skip));
  pathname <- AffymetrixFile$renameToUpperCaseExt(pathname);

  # Already corrected?
  if (isFile(pathname) && skip) {
    verbose && cat(verbose, "Optical background adjusted data file already exists: ", pathname);
    # make sure CDF gets inherited from source object
    res <- fromFile(this, pathname);
    setCdf(res, cdf);
    return(res);
  }

  # Get all probe signals
  verbose && enter(verbose, "Reading probe intensities");
  x <- getData(this, fields="intensities", verbose=less(verbose,2));
  x <- x$intensities;
  verbose && exit(verbose);

  # Identify the subset of probes to be updated?
  if (getFraction || !is.null(typesToUpdate)) {
    verbose && enter(verbose, "Identifying probes to be updated");
    subsetToUpdate <- identifyCells(cdf, indices=subsetToUpdate,
                                types=typesToUpdate, verbose=less(verbose));
    verbose && exit(verbose);
  }

  # Subtract optical background from selected probes
  verbose && enter(verbose, "Adjusting background for optical effect");
  arrayMinimum <- min(x[subsetToUpdate], na.rm=TRUE);
  verbose && printf(verbose, "Array minimum: %.2f\n", arrayMinimum);
  xdiff <- (arrayMinimum - minimum);
  verbose && printf(verbose, "Correction: -(%.2f-%.2f) = %+.2f\n", 
                                         arrayMinimum, minimum, -xdiff);
  x[subsetToUpdate] <- x[subsetToUpdate] - xdiff;
  rm(subsetToUpdate);
  verbose && exit(verbose);

  # Write adjusted data to file
  verbose && enter(verbose, "Writing adjusted probe signals");

  # Create CEL file to store results, if missing
  verbose && enter(verbose, "Creating CEL file for results, if missing");
  createFrom(this, filename=pathname, path=NULL, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing adjusted intensities");
  updateCel(pathname, intensities=x);
  verbose && exit(verbose);
  verbose && exit(verbose);

  # Return new BG adjusted data file object, making sure CDF is
  # inherited from unadjusted data
  res <- fromFile(this, pathname);
  setCdf(res, cdf);

  res;
}, private=TRUE)


###########################################################################/**
# @RdocMethod bgAdjustGcrma
#
# @title "Applies probe sequence based background correction to a CEL file"
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
#   \item{type}{The type of background correction.  Currently accepted types
#       are "fullmodel" (the default, uses MMs) and "affinities" (uses
#       probe sequence only).}
#   \item{indicesNegativeControl}{Locations of any negative control
#       probes (e.g., the anti-genomic controls on the human exon array).
#       If @NULL and type=="affinities", MMs are used as the negative
#       controls.}
#   \item{affinities}{A @numeric @vector of probe affinities, usually as
#       calculated by \code{computeAffinities()} of the 
#       @see "AffymetrixCdfFile" class.}
#   \item{gsbAdjust}{Should we adjust for specific binding (defaults to
#        @TRUE)?}
#   \item{parametersGsb}{Specific binding parameters as estimated by
#       \code{calculateParametersGsb()} for the @see "AffymetrixCelSet"
#       class.}
#   \item{k}{Tuning parameter passed to \code{bg.adjust.gcrma}.}
#   \item{rho}{Tuning parameter passed to \code{bg.adjust.gcrma}.}
#   \item{stretch}{Tuning parameter passed to \code{bg.adjust.gcrma}.}
#   \item{fast}{If @TRUE, an ad hoc transformation of the PM is performed
#       (\code{gcrma.bg.transformation.fast()}).}
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
#   Ken Simpson (ksimpson[at]wehi.edu.au) and Mark Robinson.
# }
#
# \seealso{
#  @see "gcrma::bg.adjust.gcrma"
#  @seeclass
# }
#*/###########################################################################
setMethodS3("bgAdjustGcrma", "AffymetrixCelFile", function(this, path=NULL, type="fullmodel", indicesNegativeControl=NULL, affinities=NULL, gsbAdjust=TRUE, parametersGsb=NULL, k=6*fast + 0.5*(1-fast), rho=0.7, stretch=1.15*fast + 1*(1-fast), fast=TRUE, overwrite=FALSE, skip=!overwrite, ..., verbose=FALSE, .deprecated=TRUE) {
  if (.deprecated) {
    throw("bgAdjustGcrma() is deprecated.  Please use the GcRmaBackgroundCorrection class");
  }

  # require("gcrma") || throw("Package not loaded: gcrma");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePathname(path);
  if (identical(getPath(this), path)) {
    throw("Cannot background correct data file. Argument 'path' refers to the same path as the path of the data file to be normalized: ", path);
  }

  
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- basename(getPathname(this));
  filename <- gsub("[.]cel$", ".CEL", filename);  # Only output upper case!
  pathname <- Arguments$getWritablePathname(filename, path=path,
                                            mustNotExist=(!overwrite && !skip));
  
  cdf <- getCdf(this);

  # Already corrected?
  if (isFile(pathname) && skip) {
    verbose && cat(verbose, "GC-adjusted data file already exists: ", pathname);
    # inheritance of CDF
    res <- fromFile(this, pathname);
    setCdf(res, cdf);
    return(res);
  }
  
  if (is.null(affinities)) {
# try to find APD file containing probe affinities

    paths <- path;
    paths <- paste(".",
                   paths,
                   "data/",
                   sep=";", collapse=";");

    pattern <- paste(getChipType(getCdf(this)), "-affinities.apa", sep="");
    affinityFilename <- findFiles(pattern=pattern, paths=paths, firstOnly=TRUE);
    if (is.null(affinityFilename))
      throw("Could not locate probe affinities file: ", pattern);
    affinities <- readApd(affinityFilename)$affinities;
  }

  verbose && enter(verbose, "Obtaining PM and MM signals and affinities");
  
#  mmi <- identifyCells(cdf, types="mm");
#  pmi <- identifyCells(cdf, types="pm");
#  this should work, but gives inconsistent results with bg.adjust.gcrma().
#  Stick with version below for now until we work out what is causing
#  the inconsistency.
  
  chipType <- getChipType(this);
  key <- list(method="bgAdjustGcrma", class=class(this)[1], chipType=chipType, source="gcrma");
  dirs <- c("aroma.affymetrix", chipType);
  indices <- loadCache(key=key, dirs=dirs);
  if (is.null(indices)) {
    indices <- getCellIndices(cdf, useNames=FALSE, unlist=TRUE);
  }
  saveCache(indices, key=key, dirs=dirs);
  
  # PM and MM
  mm <- getData(this, indices=indices[!isPm(cdf)])$intensities;
  pm <- getData(this, indices=indices[isPm(cdf)])$intensities;

  # corresponding affinities
  #apm <- affinities[indices[!isPm(cdf)]];
  #amm <- affinities[indices[isPm(cdf)]];
  apm <- affinities[indices[isPm(cdf)]];
  amm <- affinities[indices[!isPm(cdf)]];

  verbose && exit(verbose);

  verbose && enter(verbose, "Estimating non-specific binding parameters")

  # affinity and intensity for negative control probes
  anc <- NULL;
  ncs <- NULL;

  pmKeep <- NULL;
  
  if (!is.null(indicesNegativeControl)) {
    anc <- affinities[indicesNegativeControl]
    ncs <- getData(this, indices=indicesNegativeControl)$intensities;
  }
  
  # adjust background - use original GCRMA functions to avoid errors from
  # re-coding
  if (type == "fullmodel") {
    pm <- gcrma::bg.adjust.fullmodel(pm, mm, ncs=ncs, apm, amm, anc=anc, index.affinities=1:length(pm), k=k, rho=rho, fast=fast);
  } else if (type == "affinities") {
    if (is.null(ncs)) {
      # use MM as negative controls
      pm <- gcrma::bg.adjust.affinities(pm, mm, apm, amm, index.affinities=1:length(pm), k=k, fast=fast);
    } else {
      # use specified negative controls
      keep <- whichVector(!is.na(anc) & !is.na(ncs));
      anc <- anc[keep];
      ncs <- ncs[keep];
      pmKeep <- whichVector(!is.na(apm) & !is.na(pm));
      apm <- apm[pmKeep];
      pm <- pm[pmKeep];
      pm <- gcrma::bg.adjust.affinities(pm, ncs, apm, anc, index.affinities=1:length(pm), k=k, fast=fast, nomm=TRUE);
    }
  }
    
  verbose && exit(verbose);

  # if specific binding correction requested, carry it out
  if (gsbAdjust && !is.null(parametersGsb)) {
    verbose && enter(verbose, "Adjusting for specific binding")
	         #> GSB.adj
	         #function (Yin, subset, aff, fit1, k = k) 
	         #{
	         #    y0 = Yin[subset]
	         #    y0.adj = k + 2^(-fit1[2] * (aff - mean(aff))) * (y0 - k)
	         #    Yin[subset] = y0.adj
	         #    Yin
	         #}
    #pm <- 2^(log2(pm) - parametersGsb[2]*apm + mean(parametersGsb[2]*apm));  # this is what it used to be
	pm <- k + 2^(-parametersGsb[2] * (apm - mean(apm,na.rm=TRUE))) * (pm - k)
    verbose && exit(verbose);
  }

  # don't understand this, but it was in original bg.adjust.gcrma, so
  # we will keep it
  if (stretch != 1) {
    mu <- mean(log(pm), na.rm=TRUE);
    pm <- exp(mu + stretch * (log(pm) - mu));
  }
    
  
  # update the PM

  # Write adjusted data to file
  verbose && enter(verbose, "Writing adjusted probe signals");

  # Create CEL file to store results, if missing
  verbose && enter(verbose, "Creating CEL file for results, if missing");
  createFrom(this, filename=pathname, path=NULL, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing adjusted intensities");
  cells <- indices[isPm(cdf)];
  if (!is.null(pmKeep)) {
    cells <- cells[pmKeep];
  }
  updateCel(pathname, indices=cells, intensities=pm);
  verbose && exit(verbose);
  verbose && exit(verbose);

  # Return new background corrected data file object

  # inheritance of CDF
  res <- fromFile(this, pathname);
  setCdf(res, cdf);

  res;
}, private=TRUE)


###########################################################################/**
# @RdocMethod bgAdjustRma
#
# @title "Applies normExp background correction to a CEL file"
#
# \description{
#  @get "title".
#
#  Calls @see "affy::bg.adjust" from the \pkg{affy} package.
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The path where to save the adjusted data files.}
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
setMethodS3("bgAdjustRma", "AffymetrixCelFile", function(this, path=NULL, pmonly=TRUE, addJitter=FALSE, jitterSd=0.2, overwrite=FALSE, skip=!overwrite, ..., verbose=FALSE, .deprecated=TRUE) {
  if (.deprecated) {
    throw("bgAdjustRma() is deprecated.  Please use the RmaBackgroundCorrection class");
  }

  # Load required packages
  require("affy") || throw("Package not loaded: affy");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePathname(path);
  if (identical(getPath(this), path)) {
    throw("Cannot background correct data file. Argument 'path' refers to the same path as the path of the data file to be normalized: ", path);
  }

  
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- basename(getPathname(this));
  filename <- gsub("[.]cel$", ".CEL", filename);  # Only output upper case!
  pathname <- Arguments$getWritablePathname(filename, path=path,
                                            mustNotExist=(!overwrite && !skip));
  
  cdf <- getCdf(this);
  setCdf(this, cdf);
  
  # Already corrected?
  if (isFile(pathname) && skip) {
    verbose && cat(verbose, "Background-adjusted data file already exists: ", pathname);
    # inheritance of CDF
    res <- fromFile(this, pathname);
    setCdf(res, cdf);
    return(res);
  }

  if (pmonly) {
    verbose && cat(verbose, "Adjusting PM signals only");
  }
  
  verbose && enter(verbose, "Obtaining signals");
  
  if (pmonly) {
    chipType <- getChipType(this);
    key <- list(method="bgAdjustRma", class=class(this)[1], chipType=chipType);
    dirs <- c("aroma.affymetrix", chipType);
    pmi <- loadCache(key=key, dirs=dirs);
    if (is.null(pmi)) {
      indices <- getCellIndices(cdf, useNames=FALSE, unlist=TRUE);
      pmi <- indices[isPm(cdf)];
    }
    saveCache(pmi, key=key, dirs=dirs);
  } else {
    pmi <- NULL;
  }
  
  pm <- getData(this, indices=pmi, "intensities")$intensities;
  if (addJitter) {
    set.seed(6022007);
    pm <- pm + rnorm(length(pm), mean=0, sd=jitterSd);
  }
  clearCache(this);
  
  verbose && exit(verbose);

  # adjust background - use original affy functions to avoid errors from
  # re-coding
  verbose && enter(verbose, "Applying normal+exponential signal model");
  pm <- bg.adjust(pm);  # From package 'affy' (without a namespace)
  verbose && exit(verbose);
  
  # update the PM

  # Write adjusted data to file
  verbose && enter(verbose, "Writing adjusted probe signals");

  # Create CEL file to store results, if missing
  verbose && enter(verbose, "Creating CEL file for results, if missing");
  createFrom(this, filename=pathname, path=NULL, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing adjusted intensities");
  updateCel(pathname, indices=pmi, intensities=pm);
  verbose && exit(verbose);

  verbose && exit(verbose);

  # get rid of redundant objects to save space
  rm(pm); rm(pmi);

  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

  # Return new background corrected data file object

  # inheritance of CDF
  res <- fromFile(this, pathname);
  setCdf(res, cdf);

  res;
}, private=TRUE)


############################################################################
# HISTORY:
# 2009-04-09 [MR]
# BUG FIX: fixed discrepancy b/w aroma.affymetrix's gene specific binding
# adjustment and gcrma's GSB.adj
# 2009-03-29 [MR]
# o Made slight modifications for bgAdjustGcRma() to work with the 
#   newer Gene 1.0 ST arrays.
# 2007-12-08
# o Now bgAdjustRma() no longer assumes 'affy' is loaded.
# 2007-06-30
# o Added .deprecated=TRUE to all methods.
# 2007-03-26
# o Added verbose output about estimated background.
# 2007-03-23
# o Replaced all usage of copyCel() with createFrom().  This allow us to 
#   later update createFrom() to create CEL files of different versions.
# 2007-03-22
# o rename gsbParameters to parametersGsb to avoid clash of arguments
#   in bgAdjustGcrma.AffymetrixCelFile().  Not sure why gsbAdjust and
#   gsbParameters are being matched, but there you go.
# 2006-10-10
# o add RMA background correction
# 2006-10-06
# o make sure cdf association is inherited
# 2006-10-04
# o Debugged, tested for consistency with bg.adjust.gcrma(), docs added
# 2006-09-28
# o Created (based on AffymetrixCelFile.normalizeQuantile.R)
############################################################################
