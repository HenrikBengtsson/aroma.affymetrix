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
# \author{Ken Simpson (ksimpson[at]wehi.edu.au).}
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
  params <- NextMethod(generic="getParameters", object=this, ...);

  pmOnly <- (this$.typesToUpdate=="pm");
  
  # Get parameters of this class
  params2 <- list(
    addJitter = this$.addJitter,
    jitterSd = this$.jitterSd,
    pmonly = pmOnly
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, private=TRUE)



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
#  Returns a @double @vector.
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
    return(invisible(outputDataSet));
  }
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get the output path
  outputPath <- getPath(this);

  args <- c(list(ds, path=outputPath, verbose=verbose, overwrite=force), params, .deprecated=FALSE);

  outputDataSet <- do.call("bgAdjustRma", args=args);
  
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  # Update the output data set
  this$.outputDataSet <- outputDataSet;

  outputDataSet;
})



############################################################################
# HISTORY:
# 2007-06-30
# o Added Rdoc comments about jitter.
# 2007-05-26
# o Updated the Rdocs.
# 2007-03-21
# o Created.
############################################################################
