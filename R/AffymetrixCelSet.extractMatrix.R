###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod extractMatrix
#
# @title "Extract data as a matrix for a set of arrays"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{cells}{(The subset of cells to be matched.
#     If @NULL, all cells are considered.}
#   \item{...}{Not used.}
#   \item{field}{The field to be extracted.}
#   \item{drop}{If @TRUE, singleton dimensions are dropped.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an JxK @double @matrix where J is the number of units,
#  and K is the number of arrays.
#  The names of the columns are the names of the arrays.
#  No names are set for the rows.
#  The rows are ordered according to \code{cells} argument.
# }
#
# @author "HB, MR"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractMatrix", "AffymetrixCelSet", function(this, cells=NULL, ..., field=c("intensities", "stdvs", "pixels"), drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  cdf <- getCdf(this)
  if (is.null(cells)) {
    ncells <- nbrOfCells(cdf)
  } else {
    cells <- Arguments$getIndices(cells, max=nbrOfCells(cdf))
    ncells <- length(cells)
  }

  # Argument 'field':
  field <- match.arg(field)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  # Settings
  gcArrayFrequency <- getOption(aromaSettings, "memory/gcArrayFrequency")
  if (is.null(gcArrayFrequency))
    gcArrayFrequency <- 10


  verbose && enter(verbose, "Getting data for the array set")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate return matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocating matrix")
  arrayNames <- getNames(this)
  nbrOfArrays <- length(arrayNames)
  if (field %in% c("pixels")) {
    naValue <- as.integer(NA)
  } else {
    naValue <- as.double(NA)
  }
  df <- matrix(naValue, nrow=ncells, ncol=nbrOfArrays)
  colnames(df) <- arrayNames
  verbose && str(verbose, df)
  verbose && printf(verbose, "RAM: %s\n", hsize(object.size(df), digits = 2L, standard = "IEC"))
  verbose && exit(verbose)

  if (!is.null(cells)) {
    verbose && enter(verbose, "Optimize reading order")
    srt <- sort(cells, method="quick", index.return=TRUE)
    o <- srt$ix
    cells <- srt$x
    # Not needed anymore
    srt <- NULL
    verbose && exit(verbose)
  } else {
    o <- seq_len(ncells)
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving data")
  for (aa in seq_len(nbrOfArrays)) {
    verbose && printf(verbose, "Array %d,\n", aa)
    cf <- this[[aa]]
    df[o,aa] <- getData(cf, indices=cells, fields=field,
                                           verbose=less(verbose))[[field]]
    if (aa %% gcArrayFrequency == 0) {
      # Garbage collect
      gc <- gc()
      verbose && print(verbose, gc)
    }
  } # for (aa in ...)
  verbose && exit(verbose)

  # Drop singleton dimensions?
  if (drop) {
    df <- drop(df)
  }

  verbose && exit(verbose)

  df
}) # extractMatrix()
