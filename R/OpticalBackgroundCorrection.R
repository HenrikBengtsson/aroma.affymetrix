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
  params <- NextMethod("getParameters");

  # Get parameters of this class
  params2 <- list(
    minimum = this$.minimum
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
  subsetToUpdate <- params$subsetToUpdate;
  typesToUpdate <- params$typesToUpdate;
  minimum <- params$minimum;
  rm(params);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying the cells to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying cells to be updated");
  subsetToUpdate <- identifyCells(cdf, indices=subsetToUpdate,
                                                      types=typesToUpdate);
  verbose && cat(verbose, "Number of cells: ", length(subsetToUpdate));
  verbose && str(verbose, subsetToUpdate);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Background correct
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(ds);
  verbose && cat(verbose, "Number of arrays: ", nbrOfArrays);
  for (ii in seq_along(ds)) {
    df <- getFile(ds, ii);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", ii, getName(df), nbrOfArrays));

    dfD <- bgAdjustOptical(df, path=outputPath, subsetToUpdate=subsetToUpdate, typesToUpdate=NULL, minimum=minimum, overwrite=force, verbose=less(verbose), .deprecated=FALSE);
    verbose && print(verbose, dfD);

    # Not needed anymore
    rm(df, dfD);

    verbose && exit(verbose);
  } # for (ii ...)

  rm(subsetToUpdate);

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
# o CLEANUP: process() for OpticalBackgroundCorrection now processes
#   each file by itself, i.e. it no longer calls bgAdjustOptical() for
#   AffymetrixCelSet (which has been removed).
# 2007-08-24
# o BUG FIX: Forgot to pass argument '.deprecated=FALSE' to bgAdjustGcrma()
#   because the latter is deprecated at the user-API level.
# 2007-03-22
# o Created.
############################################################################
