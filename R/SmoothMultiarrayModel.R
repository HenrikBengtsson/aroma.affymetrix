###########################################################################/**
# @RdocClass SmoothMultiarrayModel
#
# @title "The SmoothMultiarrayModel class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a chromosomal smoothing method done
#  chromosome by chromosome.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of @see "ChromosomalModel".}
#   \item{typoOfWeights}{A @character string.}
#   \item{bandwidth}{A single @numeric specifying the smoothing bandwidth 
#     in units of nucleotides.}
#   \item{tags}{A @character @vector of tags to be added.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#  @see "CopyNumberSegmentationModel".
# }
#*/###########################################################################
setConstructorS3("SmoothMultiarrayModel", function(..., typoOfWeights=c("none", "1/s2"), bandwidth=10e3, tags="*") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'typoOfWeights':
  typoOfWeights <- match.arg(typoOfWeights);

  # Argument 'bandwidth':
  bandwidth <- Arguments$getDouble(bandwidth, range=c(1,Inf));


  extend(ChromosomalModel(..., tags=tags), "SmoothMultiarrayModel",
    .outTuple = NULL,
    .shift = 0,
    .kernel = "gauss",
    .typoOfWeights = typoOfWeights,
    .bandwidth = bandwidth
  )
}, abstract=TRUE)



setMethodS3("as.character", "SmoothMultiarrayModel", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", this, ...);
  s <- c(s, sprintf("Kernel: %s", this$.kernel));
  s <- c(s, sprintf("Bandwidth: %.2fkb", getBandwidth(this)/1e3));
  class(s) <- "GenericSummary";
  s;
}, protected=TRUE)


setMethodS3("clearCache", "SmoothMultiarrayModel", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c()) {
    this[[ff]] <- NULL;
  }

  if (!is.null(this$.outTuple)) {
   this$.outTuple <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
})


setMethodS3("getAsteriskTags", "SmoothMultiarrayModel", function(this, collapse=NULL, ...) {
  classTag <- toupper(gsub("Model$", "", class(this)[1]));
  weightsTag <- switch(this$.typoOfWeights, "1/s2"="w=s2inv", "");
  kernelTag <- this$.kernel;
  bandwidthTag <- sprintf("b=%d", getBandwidth(this));
  tags <- c(classTag, weightsTag, kernelTag, bandwidthTag);
  tags <- tags[nchar(tags) > 0];

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


setMethodS3("getRootPath", "SmoothMultiarrayModel", function(this, ...) {
  tag <- getAsteriskTags(this)[1];
  sprintf("%sData", tolower(tag));
})

setMethodS3("getBandwidth", "SmoothMultiarrayModel", function(this, ...) {
  this$.bandwidth;
})

setMethodS3("setBandwidth", "SmoothMultiarrayModel", function(this, bandwidth, ...) {
  # Argument 'bandwidth':
  bandwidth <- Arguments$getDouble(bandwidth, range=c(1,Inf));

  oldValue <- this$.bandwidth;
  this$.bandwidth <- bandwidth;
  invisible(oldValue);
})

setMethodS3("getOutputTuple", "SmoothMultiarrayModel", function(this, ..., force=FALSE,  vebose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  outTuple <- this$.outTuple;
  if (force || is.null(outTuple)) {
    outTuple <- createOutputTuple(this, force=force, verbose=less(verbose, 2));
    this$.outTuple <- outTuple;
  }

  outTuple;
})


##############################################################################
# HISTORY:
# 2007-09-20
# o Created.
##############################################################################
