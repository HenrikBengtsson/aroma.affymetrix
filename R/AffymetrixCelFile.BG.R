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
setMethodS3("bgAdjustOptical", "AffymetrixCelFile", function(this, path, minimum=1, subsetToUpdate=NULL, typesToUpdate=NULL, overwrite=FALSE, skip=!overwrite, verbose=FALSE, ..., .deprecated=TRUE) {
  if (.deprecated) {
    throw("bgAdjustOptical() is deprecated.  Please use the OpticalBackgroundCorrection class");
  }

  if (is.null(path)) {
    throw("DEPRECATED: Argument 'path' to bgAdjustOptical() must no longer be NULL.");
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
    subsetToUpdate <- Arguments$getIndices(subsetToUpdate, max=nbrOfCells);
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
#       If @NULL and \code{type == "affinities"}, then MMs are used as
#       the negative controls.}
#   \item{affinities}{A @numeric @vector of probe affinities, usually as
#       calculated by \code{computeAffinities()} of the 
#       @see "AffymetrixCdfFile" class.}
#   \item{gsbAdjust}{If @TRUE, adjustment for specific binding is done,
#       otherwise not.}
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
setMethodS3("bgAdjustGcrma", "AffymetrixCelFile", function(this, path, type=c("fullmodel", "affinities"), indicesNegativeControl=NULL, affinities=NULL, gsbAdjust=TRUE, parametersGsb=NULL, k=ifelse(fast,6,0.5), rho=0.7, stretch=ifelse(fast,1.15,1), fast=TRUE, overwrite=FALSE, skip=!overwrite, ..., verbose=FALSE, .deprecated=TRUE) {
  if (.deprecated) {
    throw("bgAdjustGcrma() is deprecated.  Please use the GcRmaBackgroundCorrection class");
  }

  if (is.null(path)) {
    throw("DEPRECATED: Argument 'path' to bgAdjustGcrma() must no longer be NULL.");
  }

  if (is.null(affinities)) {
    throw("DEPRECATED: bgAdjustGcrma() must not be called with affinities=NULL.");
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

  # Argument 'type':
  type <- match.arg(type);

  # Argument 'indicesNegativeControl':
  if (!is.null(indicesNegativeControl)) {
    indicesNegativeControl <- Arguments$getIndices(indicesNegativeControl, range=c(1, nbrOfCells(this)));
  }

  # Argument 'gsbAdjust':
  gsbAdjust <- Arguments$getLogical(gsbAdjust);

  # Argument 'k':
  k <- Arguments$getNumeric(k);

  # Argument 'rho':
  rho <- Arguments$getNumeric(rho);

  # Argument 'stretch':
  stretch <- Arguments$getNumeric(stretch);

  # Argument 'fast':
  fast <- Arguments$getLogical(fast);
  
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

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading probe signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving PM and MM signals and affinities");
 
  chipType <- getChipType(this);
  key <- list(method="bgAdjustGcrma", class=class(this)[1], chipType=chipType, source="gcrma");
  dirs <- c("aroma.affymetrix", chipType);
  indices <- loadCache(key=key, dirs=dirs);
  if (is.null(indices)) {
    indices <- getCellIndices(cdf, useNames=FALSE, unlist=TRUE);
  }
  saveCache(indices, key=key, dirs=dirs);

  # Identify PM & MM cell indices
  # Ordered according to CEL file [whilst isPm() is ordered as the CDF]
  pmCells <- indices[isPm(cdf)];
  mmCells <- indices[!isPm(cdf)];
  verbose && cat(verbose, "Number of PM probes: ", length(pmCells));
  verbose && cat(verbose, "Number of MM probes: ", length(mmCells));
  
  # PM & MM signals
  pm <- getData(this, indices=pmCells)$intensities;
  verbose && cat(verbose, "PM signals:");
  verbose && str(verbose, pm);
  mm <- getData(this, indices=mmCells)$intensities;
  verbose && cat(verbose, "MM signals:");
  verbose && str(verbose, mm);

  # PM & MM affinities
  apm <- affinities[pmCells];
  amm <- affinities[mmCells];
  verbose && cat(verbose, "PM affinities:");
  verbose && str(verbose, apm);
  verbose && cat(verbose, "MM affinities:");
  verbose && str(verbose, amm);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimating non-specific binding parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Estimating non-specific binding parameters")
  verbose && cat(verbose, "Model type/flavor: ", type);

  # affinity and intensity for negative control probes
  anc <- NULL;
  ncs <- NULL;

  if (!is.null(indicesNegativeControl)) {
    verbose && enter(verbose, "Retreaving affinities and probe signals for specified negative controls")
    anc <- affinities[indicesNegativeControl];
    verbose && cat(verbose, "Probe affinities for negative controls:");
    verbose && str(verbose, anc);
    ncs <- getData(this, indices=indicesNegativeControl)$intensities;
    verbose && cat(verbose, "Probe signals for negative controls:");
    verbose && str(verbose, ncs);
    stopifnot(length(ncs) == length(anc));
    verbose && exit(verbose);
  }
  
  # adjust background - use original GCRMA functions to avoid errors from
  # re-coding
  if (type == "fullmodel") {
    verbose && enter(verbose, "Full GCRMA model background adjustment");

    verbose && cat(verbose, "Number of PMs: ", length(pm));
    verbose && cat(verbose, "Number of MMs: ", length(mm));

    pm <- gcrma::bg.adjust.fullmodel(pms=pm, mms=mm, ncs=ncs, apm=apm, amm=amm, anc=anc,
                      index.affinities=seq_len(length(pm)), k=k, rho=rho, fast=fast);

    verbose && exit(verbose);
  } else if (type == "affinities") {
    verbose && enter(verbose, "Affinity-based GCRMA model background adjustment");

    if (is.null(ncs)) {
      verbose && cat(verbose, "Using mismatch probes (MMs) as negative controls");
      ncs <- mm;
      anc <- amm;
      nomm <- FALSE;  # However...
      # AD HOC: If there are equal number of PMs and MMs, then we assume they
      # are matched, and we will allow gcrma::bg.adjust.affinities() to subset 
      # also MMs using 'index.affinities', otherwise not.  See its code:
      #  if (!nomm)
      #    parameters <- bg.parameters.ns(ncs[index.affinities], anc, apm)
      #  else
      #    parameters <- bg.parameters.ns(ncs, anc, apm)
      # However, since we use 'index.affinities' to use all PMs, this special
      # case would equal using all MMs as well, and in that case we can equally
      # well use nomm=TRUE (see its code).
      # /HB 2010-10-02
      nomm <- TRUE;
    } else {
      verbose && cat(verbose, "Using a specified set of negative controls");
      nomm <- TRUE;
    }

    verbose && enter(verbose, "Dropping perfect-match probes (PMs) with missing signals or missing affinities");
    n0 <- length(pm);
    keepA <- (!is.na(pm));
    n <- sum(keepA);
    verbose && printf(verbose, "Number of finite PMs: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0);

    keepB <- (!is.na(apm));
    n <- sum(keepB);
    verbose && printf(verbose, "Number of finite PM affinities: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0);

    keep <- which(keepA & keepB);
    n <- length(keep);
    verbose && printf(verbose, "Number of PMs kept: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0);

    pm <- pm[keep];
    pmCells <- pmCells[keep];
    apm <- apm[keep];
    verbose && exit(verbose);
    rm(keep);

    verbose && enter(verbose, "Dropping negative controls with missing signals or missing affinities");
    n0 <- length(ncs);
    keepA <- (!is.na(ncs));
    n <- sum(keepA);
    verbose && printf(verbose, "Number of finite negative controls: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0);

    keepB <- (!is.na(anc));
    n <- sum(keepB);
    verbose && printf(verbose, "Number of finite negative-control affinities: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0);

    keep <- which(keepA & keepB);
    n <- length(keep);
    verbose && printf(verbose, "Number of negative controls kept: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0);

    ncs <- ncs[keep];
    anc <- anc[keep];
    verbose && exit(verbose);
    rm(keep);

    verbose && cat(verbose, "Number of PMs: ", length(pm));
    verbose && cat(verbose, "Number of negative controls: ", length(ncs));

    minNbrOfNegControls <- 1L;
    if (length(ncs) < minNbrOfNegControls) {
      throw(sprintf("Cannot perform GCRMA background (type=\"affinities\") correction: The number (%d) of negative control is too small.", length(ncs)));
    }

    pm <- gcrma::bg.adjust.affinities(pms=pm, ncs=ncs, apm=apm, anc=anc, 
          index.affinities=seq_len(length(pm)), k=k, fast=fast, nomm=nomm);

    verbose && exit(verbose);
  } # if (type == ...)

  # Not needed anymore
  rm(anc, ncs, mm, amm);

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
    pm <- k + 2^(-parametersGsb[2] * (apm - mean(apm, na.rm=TRUE))) * (pm - k);
    verbose && exit(verbose);
  }

  # Not needed anymore
  rm(apm);


  # don't understand this, but it was in original bg.adjust.gcrma(), so
  # we will keep it. /KS
  if (stretch != 1) {
    verbose && enter(verbose, "Stretching");
    verbose && cat(verbose, "Stretch factor: ", stretch);
    mu <- mean(log(pm), na.rm=TRUE);
    pm <- exp(mu + stretch * (log(pm) - mu));
    verbose && exit(verbose);
  }
    
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store the adjusted PM signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write adjusted data to file
  verbose && enter(verbose, "Writing adjusted probe signals");

  # Create CEL file to store results, if missing
  verbose && enter(verbose, "Creating CEL file for results, if missing");
  createFrom(this, filename=pathname, path=NULL, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing adjusted intensities");
  verbose && cat(verbose, "Number of cells (PMs only): ", length(pmCells));
  updateCel(pathname, indices=pmCells, intensities=pm);
  verbose && exit(verbose);
  verbose && exit(verbose);

  # Return new background corrected data file object

  # inheritance of CDF
  res <- fromFile(this, pathname);
  setCdf(res, cdf);

  res;
}, private=TRUE)  # bgAdjustGcrma()


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
setMethodS3("bgAdjustRma", "AffymetrixCelFile", function(this, path, pmonly=TRUE, addJitter=FALSE, jitterSd=0.2, overwrite=FALSE, skip=!overwrite, ..., verbose=FALSE, .deprecated=TRUE) {
  if (.deprecated) {
    throw("bgAdjustRma() is deprecated.  Please use the RmaBackgroundCorrection class");
  }

  if (is.null(path)) {
    throw("DEPRECATED: Argument 'path' to bgAdjustRma() must no longer be NULL.");
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
    pmCells <- loadCache(key=key, dirs=dirs);
    if (is.null(pmCells)) {
      indices <- getCellIndices(cdf, useNames=FALSE, unlist=TRUE);
      pmCells <- indices[isPm(cdf)];
    }
    saveCache(pmCells, key=key, dirs=dirs);
  } else {
    pmCells <- NULL;
  }
  
  pm <- getData(this, indices=pmCells, "intensities")$intensities;
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
  updateCel(pathname, indices=pmCells, intensities=pm);
  verbose && exit(verbose);

  verbose && exit(verbose);

  # get rid of redundant objects to save space
  rm(pm, pmCells);

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
# 2012-11-20
# o CLEANUP: Removed obsolete code from internal bgAdjustOptical() that 
#   was never reached and that loaded affinities via obsolete APD files.
# 2011-03-26
# o ROBUSTNESS: Now internal bgAdjustGcrma(..., type="affinities") for
#   AffymetrixCelFile gives a more informative error message when there
#   are too few negative controls.
# 2010-10-02
# o We now use nomm=TRUE for all cases for type="affinities".  See code
#   for explanation.  This also solves the problems of using for instance 
#   chip type MoEx-1_0-st-v1.
# 2010-09-29
# o ROBUSTNESS: Now bgAdjustGcrma(..., affinities=NULL) is deprecated and
#   throws an exception.
# o CLEANUP: Cleaned up bgAdjustGcrma().
# o Added more verbose output to bgAdjustGcrma() for AffymetrixCelFile.
# 2009-04-09 [MR]
# o BUG FIX: fixed discrepancy b/w aroma.affymetrix's gene specific binding
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
