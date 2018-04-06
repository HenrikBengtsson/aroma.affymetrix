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
# @author "HB"
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
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Getting CEL file to store baseline signals")
  key <- list(dataset=getFullName(this), samples=getNames(this))
  id <- getChecksum(key)
  filename <- sprintf(".baseline,%s.CEL", id)


  verbose && enter(verbose, "Searching for an existing file")

  # Searching for the output file in multiple directories
  path <- getPath(this)
  paths <- c(path)

  # Drop tags from root path?
  if (getOption(aromaSettings, "devel/dropRootPathTags", TRUE)) {
    path <- dropRootPathTags(path, depth=2, verbose=less(verbose, 5))
    paths <- c(paths, path)
    paths <- unique(paths)
  }

  verbose && cat(verbose, "Paths:")
  verbose && print(verbose, paths)
  verbose && cat(verbose, "Filename: ", filename)

  pathname <- NULL
  for (kk in seq_along(paths)) {
    path <- paths[kk]
    verbose && enter(verbose, sprintf("Searching path #%d of %d", kk, length(paths)))

    verbose && cat(verbose, "Path: ", path)
    pathnameT <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE)
    verbose && cat(verbose, "Pathname: ", pathnameT)
    if (isFile(pathnameT)) {
      pathname <- pathnameT
      verbose && cat(verbose, "Found an existing file.")
      verbose && exit(verbose)
      break
    }

    verbose && exit(verbose)
  } # for (kk ...)
  verbose && cat(verbose, "Located pathname: ", pathname)

  verbose && exit(verbose)


  # Get a template CEL file
  df <- getOneFile(this)

  if (isFile(pathname)) {
    verbose && enter(verbose, "Loading existing data file")
    res <- newInstance(df, pathname)
    verbose && exit(verbose)
  } else {
    verbose && enter(verbose, "Allocating empty data file")

    path <- paths[length(paths)]
    verbose && cat(verbose, "Path: ", path)
    verbose && cat(verbose, "Filename: ", filename)
    pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=TRUE)

    verbose && enter(verbose, "Retrieving CEL file")
    res <- createFrom(df, filename=pathname, path=NULL, methods="create",
                         clear=TRUE, force=force, verbose=less(verbose))
    verbose && print(verbose, res)
    verbose && exit(verbose)

    verbose && exit(verbose)
  } # if (isFile(pathname))

  # Not needed anymore
  df <- NULL
  verbose && exit(verbose)

  res
}, protected=TRUE)



setMethodS3("getBaseline", "SnpChipEffectSet", function(this, ...) {
  res <- NextMethod("getBaseline")
  res$mergeStrands <- getMergeStrands(this)
  res
}, protected=TRUE)


setMethodS3("getBaseline", "CnChipEffectSet", function(this, ...) {
  res <- NextMethod("getBaseline")
  res$combineAlleles <- getCombineAlleles(this)
  res
}, protected=TRUE)
