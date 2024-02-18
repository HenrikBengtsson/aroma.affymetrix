setConstructorS3("Transform", function(..., .reqSetClass="AffymetrixCelSet") {
  extend(AromaTransform(..., .reqSetClass=.reqSetClass), "Transform")
}, abstract=TRUE)


setMethodS3("getOutputFiles", "Transform", function(this, pattern=NULL, ...) {
  # Argument 'pattern':
  if (is.null(pattern)) {
    # Default filename pattern find non-private (no dot prefix) CEL files.
    pattern <- "^[^.].*[.](cel|CEL)$"
  } else {
    pattern <- Arguments$getRegularExpression(pattern=pattern)
  }

  NextMethod("getOutputFiles", pattern=pattern)
}, protected=TRUE)



###########################################################################/**
# @set "class=Transform"
# @RdocMethod getOutputDataSet
#
# @title "Gets the transformed data set"
#
# \description{
#  @get "title", if processed.
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
#  Returns an @see "aroma.core::AromaMicroarrayDataSet" or @NULL.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getOutputDataSet", "Transform", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Getting output data set for ", class(this)[1])

  # Inherit the CDF from the input data set.
  ds <- getInputDataSet(this)
  cdf <- getCdf(ds)
  args <- list(generic="getOutputDataSet", this, ...,
               cdf=cdf, checkChipType=FALSE)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Inherit certain arguments from the input data set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # AD HOC (not using OO), but setting these arguments does speed
  # up things. /HB 2007-09-17
  # Note, this is also done in Transform for now, such it is not really
  # needed here.  However, in case it will be removed from there it still
  # makes sense to have it here.
  if (inherits(ds, "CnChipEffectSet"))
    args$combineAlleles <- ds$combineAlleles
  if (inherits(ds, "SnpChipEffectSet"))
    args$mergeStrands <- ds$mergeStrands

  verbose && cat(verbose, "Calling NextMethod:")
  verbose && str(verbose, args)
  args$verbose <- less(verbose,1)

  res <- do.call(NextMethod, args)

  # Let the set update itself
  if (!is.null(res)) {
    update2(res, ..., verbose=less(verbose,1))
  }

  verbose && exit(verbose)

  res
})
