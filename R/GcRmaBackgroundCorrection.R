###########################################################################/**
# @RdocClass GcRmaBackgroundCorrection
#
# @title "The GcRmaBackgroundCorrection class"
#
# \description{
#  @classhierarchy
#
#  This class represents the GCRMA background adjustment function.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#      @see "ProbeLevelTransform".}
#   \item{indicesNegativeControl}{Locations of any negative control
#       probes (e.g., the anti-genomic controls on the human exon array).
#       If @NULL and \code{type == "affinities"}, then all non-PM probes
#       are used as the negative controls.}
#   \item{affinities}{A @numeric @vector of probe affinities, usually as
#       calculated by \code{computeAffinities()} of the 
#       @see "AffymetrixCdfFile" class.}
#   \item{type}{Type (flavor) of background correction, which can
#       be either \code{"fullmodel"} (uses MMs; requires that the chip type
#       has PM/MM pairs) or \code{"affinities"} (uses probe sequence only).}
#   \item{gsbAdjust}{If @TRUE, adjustment for specific binding is done,
#       otherwise not.}
#   \item{opticalAdjust}{If @TRUE, correction for optical effect is done
#       first, utilizing @see "OpticalBackgroundCorrection".}
#   \item{gsbParameters}{Additional argument passed to the internal
#       \code{bgAdjustGcrma()} method.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \author{Ken Simpson (ksimpson[at]wehi.edu.au).}
#*/###########################################################################
setConstructorS3("GcRmaBackgroundCorrection", function(..., indicesNegativeControl=NULL, affinities=NULL, type=c("fullmodel", "affinities"), opticalAdjust=TRUE, gsbAdjust=TRUE, gsbParameters=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'indicesNegativeControl':
  if (!is.null(indicesNegativeControl)) {
    indicesNegativeControl <- Arguments$getIndices(indicesNegativeControl);
  }

  # Argument 'affinities':
  if (!is.null(affinities)) {
    affinities <- Arguments$getNumerics(affinities);
  }

  # Argument 'type':
  type <- match.arg(type);

  # Argument 'opticalAdjust':
  opticalAdjust <- Arguments$getLogical(opticalAdjust);

  # Argument 'gsbAdjust':
  gsbAdjust <- Arguments$getLogical(gsbAdjust);

  # Argument 'gsbParameters':


  extend(BackgroundCorrection(..., typesToUpdate="pm"), "GcRmaBackgroundCorrection",
    .indicesNegativeControl=indicesNegativeControl,
    .affinities=affinities,
    .type=type,
    .opticalAdjust=opticalAdjust,
    .gsbAdjust=gsbAdjust,
    .gsbParameters=gsbParameters
  );
})


setMethodS3("getParameters", "GcRmaBackgroundCorrection", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  # Get parameters of this class
  params2 <- list(
    indicesNegativeControl = this$.indicesNegativeControl,
    affinities = this$.affinities,
    type = this$.type,
    opticalAdjust = this$.opticalAdjust,
    gsbAdjust = this$.gsbAdjust,
    gsbParameters = this$.gsbParameters
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
setMethodS3("process", "GcRmaBackgroundCorrection", function(this, ..., force=FALSE, verbose=FALSE) {

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

  outputDataSet <- do.call("bgAdjustGcrma", args=args);
  
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
# 2010-10-01
# o Now GcRmaBackgroundCorrection tries to calculate probe affinites based
#   on ACS annotation files and then as a backup/backward compatibility
#   it uses Affymetrix probe-tab files.
# 2010-09-26
# o Added explicit descriptions to the arguments list of the Rdocs.
# o ROBUSTNESS: Added more validation of the arguments passed to 
#   the GcRmaBackgroundCorrection constructor.
# 2007-08-24
# o BUG FIX: Forgot to pass argument '.deprecated=FALSE' to bgAdjustGcrma()
#   because the latter is deprecated at the user-API level.
# 2007-03-21
# o Created.
############################################################################

