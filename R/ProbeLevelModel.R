###########################################################################/**
# @RdocClass ProbeLevelModel
#
# @title "The ProbeLevelModel class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a probe-level model (PLM) as defined
#  by the \pkg{affyPLM} package:
#    "A [...] PLM is a model that is fit to probe-intensity data.
#     More specifically, it is where we fit a model with probe level
#     and chip level parameters on a probeset by probeset basis",
#  where the more general case for a probeset is a \emph{unit group}
#  in Affymetrix CDF terms.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "MultiArrayUnitModel".}
#   \item{standardize}{If @TRUE, chip-effect and probe-affinity estimates are
#      rescaled such that the product of the probe affinities is one.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   In order to minimize the risk for mistakes, but also to be able compare
#   results from different PLMs, all PLM subclasses must meet the following
#   criteria:
#   \enumerate{
#     \item All parameter estimates must be (stored and returned) on the
#       intensity scale, e.g. log-additive models such as @see "RmaPlm"
#       have to transform the parameters on the log-scale to the intensity
#       scale.
#     \item The probe-affinity estimates \eqn{\phi_k} for a unit group
#       must be constrained such that \eqn{\prod_k \phi_k = 1},
#       or equivalently if \eqn{\phi_k > 0},\eqn{\sum_k \log(\phi_k) = 0}.
#   }
#   Note that the above probe-affinity constraint guarantees that the
#   estimated chip effects across models are on the same scale.
# }
#
# @author "HB"
#
# \seealso{
#   For more details on probe-level models, please see
#   the \pkg{preprocessCore} package.
# }
#*/###########################################################################
setConstructorS3("ProbeLevelModel", function(..., standardize=TRUE) {
  extend(MultiArrayUnitModel(...), "ProbeLevelModel",
    "cached:.paf" = NULL,
    "cached:.ces" = NULL,
    "cached:.rs" = NULL,
    "cached:.ws" = NULL,
    "cached:.lastPlotData" = NULL,
    standardize = standardize
  )
}, abstract=TRUE)


setMethodS3("getAsteriskTags", "ProbeLevelModel", function(this, collapse=NULL, ...) {
  # Returns 'PLM' (but allow for future extensions)
  tags <- NextMethod("getAsteriskTags", collapse=NULL)
  tags[1] <- "PLM"
  tags
}, protected=TRUE)



setMethodS3("getRootPath", "ProbeLevelModel", function(this, ...) {
  "plmData"
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getProbeAffinityFile
# @aliasmethod getProbeAffinities
#
# @title "Gets the probe affinities for this model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{.class}{A @see "ProbeAffinityFile" \emph{class}.}
# }
#
# \value{
#  Returns a @see "ProbeAffinityFile" object.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getProbeAffinityFile", "ProbeLevelModel", function(this, ..., .class=ProbeAffinityFile) {
  paf <- this$.paf
  if (!is.null(paf) && isFile(getPathname(paf)))
    return(paf)

  ds <- getDataSet(this)
  if (length(ds) == 0)
    throw("Cannot create probe-affinity file. There are no CEL files in the data set.")

  # Create probe-affinity file from CEL file template
  df <- getOneFile(ds, mustExist=TRUE)
  paf <- createFrom(df, filename="probeAffinities.CEL", path=getPath(this),
                                          methods="create", clear=TRUE, ...)

  # Make it into an object of the correct class
  paf <- newInstance(.class, getPathname(paf), cdf=getCdf(ds),
                                                  probeModel=this$probeModel)

  this$.paf <- paf

  paf
})



###########################################################################/**
# @RdocMethod getChipEffectSet
# @aliasmethod getChipEffects
#
# @title "Gets the set of chip effects for this model"
#
# \description{
#  @get "title".
#  There is one chip-effect file per array.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to \code{getMonocellCdf()} of
#                                                  @see "AffymetrixCdfFile".}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "ChipEffectSet" object.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipEffectSet", "ProbeLevelModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  ces <- this$.ces
  if (!is.null(ces))
    return(ces)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create chip-effect files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Let the parameter object know about the CDF structure, because we
  # might use a modified version of the one in the CEL header.
  ds <- getDataSet(this)
  if (length(ds) == 0)
    throw("Cannot create chip-effect set. The CEL set is empty.")

  verbose && enter(verbose, "Getting chip-effect set from data set")
  # Inherit the (monocell) CDF
  cdf <- getCdf(ds)
  cdfMono <- getMonocellCdf(cdf, ..., verbose=less(verbose))

  # Gets the ChipEffects Class object
  clazz <- getChipEffectSetClass(this)
  ces <- clazz$fromDataSet(dataSet=ds, path=getPath(this), cdf=cdfMono,
                                                    verbose=less(verbose))
  verbose && exit(verbose)

  # Let the set update itself
  update2(ces, verbose=less(verbose, 1))

  # Store in cache
  this$.ces <- ces

  ces
})

setMethodS3("getChipEffectSetClass", "ProbeLevelModel", function(static, ...) {
  ChipEffectSet
}, static=TRUE, private=TRUE)




setMethodS3("getResidualSet", "ProbeLevelModel", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  rs <- this$.rs
  if (!force && !is.null(rs))
    return(rs)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create residuals files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Let the parameter object know about the CDF structure, because we
  # might use a modified version of the one in the CEL header.
  ds <- getDataSet(this)
  if (length(ds) == 0)
    throw("Cannot create residuals set. The data set is empty.")

  verbose && enter(verbose, "Getting chip-effect set from data set")
  # Gets the ResidualSet Class object
  clazz <- getResidualSetClass(this)
  rs <- clazz$fromDataSet(dataSet=ds, path=getPath(this),
                                                     verbose=less(verbose))
  # Inherit the CDF from the input data set
  cdf <- getCdf(ds)
  setCdf(rs, cdf)
  verbose && exit(verbose)

  # Store in cache
  this$.rs <- rs

  rs
})

setMethodS3("getResidualSetClass", "ProbeLevelModel", function(static, ...) {
  ResidualSet
}, static=TRUE, private=TRUE)


setMethodS3("getWeightsSet", "ProbeLevelModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  ws <- this$.ws
  if (!is.null(ws))
    return(ws)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create weights files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Let the parameter object know about the CDF structure, because we
  # might use a modified version of the one in the CEL header.
  ds <- getDataSet(this)
  if (length(ds) == 0)
    throw("Cannot create weights set. The data set is empty.")

  verbose && enter(verbose, "Getting chip-effect set from data set")
  # Gets the WeightsSet Class object
  clazz <- getWeightsSetClass(this)
  ws <- clazz$fromDataSet(dataSet=ds, path=getPath(this),
                                                     verbose=less(verbose))
  # make sure CDF is inherited
  setCdf(ws, getCdf(ds))
  verbose && exit(verbose)

  # Store in cache
  this$.ws <- ws

  ws
})

setMethodS3("getWeightsSetClass", "ProbeLevelModel", function(static, ...) {
  WeightsSet
}, static=TRUE, private=TRUE)



###########################################################################/**
# @RdocMethod findUnitsTodo
#
# @title "Identifies non-fitted units"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{verbose}{A @logical or a @see "R.utils::Verbose".}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer @vector of unit indices.
# }
#
# \seealso{
#   Internally this methods calls the same method for the
#   @see "ChipEffectSet" class.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("findUnitsTodo", "ProbeLevelModel", function(this, verbose=FALSE, ...) {
  ces <- getChipEffectSet(this, verbose=verbose)
  findUnitsTodo(ces, verbose=verbose, ...)
}, private=TRUE)
