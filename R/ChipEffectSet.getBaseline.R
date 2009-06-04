###########################################################################/**
# @set "class=ChipEffectSet"
# @RdocMethod getBaseline
#
# @title "Gets the baseline signals across chromosomes"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{force}{If @TRUE, the CEL file that stores the is recreated.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getBaseline", "ChipEffectSet", function(this, force=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting CEL file to store baseline signals");
  path <- getPath(this);
  key <- list(dataset=getFullName(this), samples=getNames(this));
  id <- digest2(key);

  # Generate output pathname
  filename <- sprintf(".baseline,%s.CEL", id);
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Get a template CEL file
  df <- getFile(this, 1);

  verbose && enter(verbose, "Retrieving CEL file");
  res <- createFrom(df, filename=pathname, path=NULL, methods="create",
                         clear=TRUE, force=force, verbose=less(verbose));
  verbose && print(verbose, res);
  verbose && exit(verbose);

  rm(df);
  verbose && exit(verbose);

  res;
})



setMethodS3("getBaseline", "SnpChipEffectSet", function(this, ...) {
  res <- NextMethod("getBaseline", this, ...);
  res$mergeStrands <- getMergeStrands(this);
  res;
})

setMethodS3("getBaseline", "CnChipEffectSet", function(this, ...) {
  res <- NextMethod("getBaseline", this, ...);
  res$combineAlleles <- getCombineAlleles(this);
  res;
})




############################################################################
# HISTORY:
# 2007-08-09
# o getBaseLine() of CnChipEffectSet now creates CEL files with upper-case
#   filename extension "*.CEL", not "*.cel".  The reason for this is that
#   some software don't recognize lower case filename extensions :(
# 2007-03-22
# o TO DO: Estimate standard errors just like getAverage() does.
# o Added getBaseline().
# o First working version of calculateBaseline(). Method now creates a CEL
#   files to store the estimates.
# 2007-03-16
# o Created.  See ploidy4.ps paper.
############################################################################
