###########################################################################/**
# @set "class=ChipEffectSet"
# @RdocMethod extractExpressionSet
#
# @title "Extracts an in-memory ExpressionSet object"
#
# \description{
#  @get "title" from a @see "ChipEffectSet" object.
#  Note that any modifications done to the extract object will \emph{not}
#  be reflected in the original chip-effect set.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Optional argument passed to @seemethod "extractMatrix".}
#   \item{logBase}{An @integer specifying the base to use when
#    log-transforming the signals.  If @NULL, the signals are not
#    transformed, but kept as is.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "Biobase::ExpressionSet-class" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("extractExpressionSet", "ChipEffectSet", function(this, ..., logBase=2, verbose=FALSE) {
  require("Biobase") || throw("Package not loaded: Biobase");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'logBase':
  if (!is.null(logBase)) {
    logBase <- Arguments$getInteger(logBase, range=c(1,Inf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extract an ExpressionSet");  

  verbose && print(verbose, this);
  verbose && cat(verbose, "Number of arrays: ", length(this));  


  verbose && enter(verbose, "Reading data");  
  Y <- extractMatrix(this, ..., returnUgcMap=TRUE, verbose=less(verbose, 5));
  ugcMap <- attr(Y, "unitGroupCellMap");
  attr(Y, "unitGroupCellMap") <- NULL;
  verbose && str(verbose, Y);
  verbose && str(verbose, ugcMap);
  verbose && exit(verbose);

  if (!is.null(logBase)) {
    verbose && enter(verbose, "Log-transforming signals");
    verbose && cat(verbose, "Log base: ", logBase);
    Y <- log(Y, base=logBase);
    verbose && str(verbose, Y);
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Creating feature names");

  verbose && enter(verbose, "Reading (unit,group) names from CDF");
  cdf <- getCdf(this);
  verbose && print(verbose, cdf);

  ugNames <- getUnitGroupNamesFromUgcMap(cdf, ugcMap=ugcMap,
                                         verbose=less(verbose, 10));
  verbose && str(verbose, ugNames);
  verbose && exit(verbose);

  # Turn (<unit>,<group>) names into <unit>,<group> names
  names <- paste(ugNames$unitName, ugNames$groupName, sep=",");
  names <- gsub(",$", "", names);
  verbose && str(verbose, names);
  verbose && exit(verbose);

  rownames(Y) <- names;

  # Not needed anymore
  rm(names, ugNames, ugcMap, cdf);

  verbose && enter(verbose, "Allocating ExpressionSet object");
  eset <- new("ExpressionSet", exprs=Y);
  verbose && print(verbose, eset);

  rm(Y); # Not needed anymore
  verbose && exit(verbose);

  verbose && exit(verbose);

  eset;
}) # extractExpressionSet()


###########################################################################
# HISTORY:
# 2011-07-14 [HB]
# o Added argument 'logBase' to extractExpressionSet().
# 2011-07-09 [HB]
# o Added extractExpressionSet() for ChipEffectSet.
# o Created.
###########################################################################
