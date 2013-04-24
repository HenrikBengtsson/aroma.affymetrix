###########################################################################/**
# @RdocClass RmaBackgroundCorrection
#
# @title "The RmaBackgroundCorrection class"
#
# \description{
#  @classhierarchy
#
#  This class represents the RMA background adjustment function.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "BackgroundCorrection".}
#   \item{addJitter}{If @TRUE, Zero-mean gaussian noise is added to the
#     signals before being background corrected.}
#   \item{jitterSd}{Standard deviation of the jitter noise added.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Jitter noise}{
#   The fitting algorithm of the RMA background correction model may not
#   converge if there too many small and discrete signals.  To overcome
#   this problem, a small amount of noise may be added to the signals
#   before fitting the model.  This is an ad hoc solution that seems to
#   work.
#   However, add Gaussian noise may generate non-positive signals.
# }
#
# \details{
#   Internally @see "affy::bg.adjust" is used to background correct the
#   probe signals.  The default is to background correct PM signals only.
# }
#
# @author "KS, HB"
#*/###########################################################################
setConstructorS3("RmaBackgroundCorrection", function(..., addJitter=FALSE, jitterSd=0.2) {
  extend(BackgroundCorrection(..., typesToUpdate="pm"),
    "RmaBackgroundCorrection",
    .addJitter=addJitter,
    .jitterSd=jitterSd
  );
})


setMethodS3("getParameters", "RmaBackgroundCorrection", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  pmOnly <- (this$.typesToUpdate == "pm");

  # Get parameters of this class
  params2 <- list(
    addJitter = this$.addJitter,
    jitterSd = this$.jitterSd,
    pmonly = pmOnly
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod process
#
# @title "Performs background correction"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the output data set.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "RmaBackgroundCorrection", function(this, ..., force=FALSE, verbose=FALSE) {

  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Background correcting data set");

  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already background corrected");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(outputDataSet);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get the output path
  outputPath <- getPath(this);

  cdf <- getCdf(ds);

  # Get algorithm parameters (including the target distribution)
  params <- getParameters(this);
  # 'subsetToUpdate' is not used and 'typesToUpdate' are used via 'pmonly'
  pmonly <- params$pmonly;
  addJitter <- params$addJitter;
  jitterSd <- params$jitterSd;
  rm(params);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Background correct
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(ds);
  verbose && cat(verbose, "Number of arrays: ", nbrOfArrays);
  for (ii in seq_along(ds)) {
    df <- getFile(ds, ii);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", ii, getName(df), nbrOfArrays));

    dfD <- bgAdjustRma(df, path=outputPath, pmonly=pmonly, addJitter=addJitter, jitterSd=jitterSd, overwrite=force, verbose=verbose, .deprecated=FALSE);
    verbose && print(verbose, dfD);

    # Not needed anymore
    rm(df, dfD);

    verbose && exit(verbose);
  } # for (ii ...)

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Get the output data set
  outputDataSet <- getOutputDataSet(this);

  verbose && exit(verbose);

  outputDataSet;
})


############################################################################
# HISTORY:
# 2012-11-20
# o CLEANUP: process() for RmaBackgroundCorrection now processes
#   each file by itself, i.e. it no longer calls bgAdjustRma() for
#   AffymetrixCelSet (which has been removed).
# 2007-06-30
# o Added Rdoc comments about jitter.
# 2007-05-26
# o Updated the Rdocs.
# 2007-03-21
# o Created.
############################################################################
