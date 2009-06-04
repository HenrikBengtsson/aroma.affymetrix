###########################################################################/**
# @RdocClass OpticalBackgroundCorrection
#
# @title "The OpticalBackgroundCorrection class"
#
# \description{
#  @classhierarchy
#
#  This class represents "optical" background adjustment.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "ProbeLevelTransform".}
#   \item{minimum}{The minimum signal allowed after adjustment.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \author{Ken Simpson (ksimpson[at]wehi.edu.au).}
#*/###########################################################################
setConstructorS3("OpticalBackgroundCorrection", function(..., minimum=1) {
  extend(BackgroundCorrection(..., typesToUpdate="pmmm"),
    "OpticalBackgroundCorrection",
    .minimum = minimum
  )
})


setMethodS3("getParameters", "OpticalBackgroundCorrection", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  # Get parameters of this class
  params2 <- list(
    minimum = this$.minimum
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
setMethodS3("process", "OpticalBackgroundCorrection", function(this, ..., force=FALSE, verbose=FALSE) {

  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Background correcting data set");

  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already background corrected for \"optical\" effects");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get algorithm parameters (including the target distribution)
  params <- getParameters(this);

  # Get the output path
  outputPath <- getPath(this);

  args <- c(list(ds, path=outputPath, verbose=verbose, overwrite=force), params, .deprecated=FALSE);

  outputDataSet <- do.call("bgAdjustOptical", args=args);
  
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
# 2007-08-24
# o BUG FIX: Forgot to pass argument '.deprecated=FALSE' to bgAdjustGcrma()
#   because the latter is deprecated at the user-API level.
# 2007-03-22
# o Created.
############################################################################
