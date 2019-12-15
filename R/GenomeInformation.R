###########################################################################/**
# @RdocClass GenomeInformation
#
# @title "The GenomeInformation class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFile".}
#   \item{.verify}{For internal use only.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("GenomeInformation", function(..., .verify=TRUE) {
  extend(GenericDataFile(...), c("GenomeInformation",
                                               uses("FileCacheKeyInterface")),
    "cached:.data"=NULL,
    "cached:.chromosomeStats"=NULL
  )
})


setMethodS3("as.character", "GenomeInformation", function(x, ...) {
  # To please R CMD check
  this <- x

  s <- NextMethod("as.character")
  s <- c(s, sprintf("Chip type: %s", getChipType(this)))
  s
}, protected=TRUE)




###########################################################################/**
# @RdocMethod verify
#
# @title "Verifies the correctness of the underlying file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (visibly) @TRUE if the file is valid, otherwise an error is
#   thrown.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("verify", "GenomeInformation", function(this, ...) {
  TRUE
}, protected=TRUE)




###########################################################################/**
# @RdocMethod byChipType
#
# @title "Static method to define a genome information set by chip type"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "GenomeInformation" object.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("byChipType", "GenomeInformation", static=TRUE, abstract=TRUE)



setMethodS3("fromDataSet", "GenomeInformation", function(static, dataSet, ...) {
  chipType <- getChipType(dataSet)
  byChipType(static, chipType=chipType, ...)
}, static=TRUE, protected=TRUE)



setMethodS3("getUnitsOnChromosome", "GenomeInformation", function(this, ...) {
  getUnitsOnChromosomes(this, ...)
})

setMethodS3("getUnitsOnChromosomes", "GenomeInformation", function(this, chromosomes, region=NULL, ..., .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument '...':
  args <- list(...)
  if ("chromosome" %in% names(args)) {
    chromosomes <- args[["chromosome"]]
  }

  allChromosomes <- getChromosomes(this)

  # Argument 'chromosomes':
  if (.checkArgs) {
    range <- range(allChromosomes)
  } else {
    range <- c(1,999)
  }
  chromosomes <- Arguments$getChromosomes(chromosomes, range=range, ...)

  # Argument 'region':
  if (.checkArgs) {
    if (!is.null(region)) {
      region <- Arguments$getNumerics(region)
      if (length(region) != 2) {
        throw("Argument 'region' must be a numeric vector of length 2: ",
                                                           length(region))
      }
      if (region[1] > region[2]) {
        throw("Argument 'region' is not ordered: c(",
                                        paste(region, collapse=", "), ")")
      }
    }
  }


  data <- getData(this)
  keep <- (data[,"chromosome"] %in% chromosomes)
  data <- data[keep,,drop=FALSE]

  # Stratify by physical position?
  if (!is.null(region)) {
    pos <- data[,"physicalPosition"]
    keep <- (!is.na(pos) & (region[1] <= pos & pos <= region[2]))
    data <- data[keep,,drop=FALSE]
  }

  # Extract the units
  units <- rownames(data)
  units <- as.integer(units)

  # CONTRACT
  # Sanity check, cf. thread 'Error with AllelicCrosstalkCalibration
  # and process' on 2011-03-03.
  units <- Arguments$getIndices(units)

  units
})


setMethodS3("readDataFrame", "GenomeInformation", function(this, ...) {
  readTableInternal(this, ...)
}, protected=TRUE)


setMethodS3("readTableInternal", "GenomeInformation", function(this, pathname, colClasses=NULL, ..., include=NULL, exclude=NULL, fill=TRUE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)

  verbose && enter(verbose, "Reading tabular data from file")
  pathname <- getPathname(this)
  verbose && cat(verbose, "Pathname: ", pathname)

  verbose && cat(verbose, "Argument 'include': ", paste(include, collapse=", "))
  verbose && cat(verbose, "Argument 'exclude': ", paste(exclude, collapse=", "))
  exclude <- setdiff(exclude, include)
  verbose && cat(verbose, "exclude\\include: ", paste(exclude, collapse=", "))
  colClasses[names(colClasses) %in% exclude] <- "NULL"
  toRead <- names(colClasses)[colClasses != "NULL"]
  verbose && cat(verbose, "Columns to be read: ", paste(toRead, collapse=", "))

  df <- readTable(pathname, colClasses=colClasses, header=TRUE, sep="\t", fill=fill, ..., verbose=less(verbose))

  colnames(df) <- toCamelCase(colnames(df))
  verbose && exit(verbose)

  df
}, private=TRUE)




setMethodS3("nbrOfUnits", "GenomeInformation", function(this, ...) {
  data <- getData(this, fields=1)
  nrow(data)
})

setMethodS3("getDataColumns", "GenomeInformation", function(this, ...) {
  data <- getData(this)
  colnames(data)
}, private=TRUE)


###########################################################################/**
# @RdocMethod getPositions
#
# @title "Gets the physical positions for a set of units"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getData".}
#  \item{na.rm}{If @TRUE, non-defined unit indices are excluded, otherwise
#      not.}
# }
#
# \value{
#   Returns an @integer @vector.
# }
#
# \seealso{
#   @seemethod "getData".
#   @seemethod "getUnitIndices".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPositions", "GenomeInformation", function(this, ..., na.rm=FALSE) {
  df <- getData(this, fields="physicalPosition", ...)
  suppressWarnings({  # Suppress warnings about NAs
    df <- as.integer(df[,1])
  })
  if (na.rm)
    df <- df[!is.na(df)]
  df
})



setMethodS3("getChromosomes", "GenomeInformation", function(this, ..., force=FALSE) {
  chromosomes <- this$.chromosomes
  if (is.null(chromosomes) || force) {
    chromosomes <- unique(getData(this, fields="chromosome")[,1])
    chromosomes <- sort(chromosomes)

    this$.chromosomes <- chromosomes
  }
  chromosomes
})


setMethodS3("getChromosomeStats", "GenomeInformation", function(this, na.rm=TRUE, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)

  verbose && enter(verbose, "Retrieving chromosome statistics")

  stats <- this$.chromosomeStats
  if (is.null(stats) || force) {
    chromosomes <- getChromosomes(this)
    nbrOfChromosomes <- length(chromosomes)
    naValue <- NA_real_
    stats <- matrix(naValue, nrow=nbrOfChromosomes, ncol=3)
    colnames(stats) <- c("min", "max", "n")
    rownames(stats) <- chromosomes
    for (kk in seq_along(chromosomes)) {
      chr <- chromosomes[kk]
      verbose && enter(verbose, sprintf("Chromosome %d ('Chr%02d') of %d",
                                          kk, chr, nbrOfChromosomes))
      pos <- getPositions(this, chromosome=chr)
      r <- range(pos, na.rm=na.rm)
      stats[kk,1:2] <- r
      stats[kk,3] <- length(pos)
      verbose && exit(verbose)
    }
    this$.chromosomeStats <- stats
  } else {
    verbose && cat(verbose, "Found cached results")
  }

  verbose && exit(verbose)

  stats
})


###########################################################################/**
# @RdocMethod plotDensity
#
# @title "Plots the density of SNPs for a given chromosome"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chromosome}{The chromsome to be displayed.}
#  \item{...}{Additional arguments passed to @seemethod "getPositions".}
#  \item{adjust}{A bandwidth parameter for the density estimator.}
#  \item{xlab}{The label on the x-axis.  If @NULL, a default will generated.}
#  \item{main}{The title of the plot.  If @NULL, a default will generated.}
#  \item{annotate}{If @TRUE, the plot is annotated with extra information.}
# }
#
# \value{
#   Returns (invisibly) the estimated density.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotDensity", "GenomeInformation", function(this, chromosome, ..., adjust=1/20, xlab=NULL, main=NULL, annotate=TRUE) {
  # Get the positions of these units
  pos <- getPositions(this, chromosome=chromosome, ...)

  # Estimate the SNP density
  d <- density(pos/10^6, from=0, adjust=adjust)

  # Plot the density
  if (is.null(xlab))
    xlab <- "Physical position (Mb)"
  if (is.null(main))
    main <- sprintf("SNP density on chromsome %s", chromosome)
  plot(d, lwd=2, xlab=xlab, main=main)

  # Annotate?
  if (annotate) {
    stext(sprintf("%d SNPs", d$n), side=3, pos=1, line=-1, cex=0.8)
    stextChipType(chipType=getChipType(this))
  }

  invisible(d)
})



###########################################################################/**
# @RdocMethod getUnitIndices
#
# @title "Gets unit indices ordered along the chromosome"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getData".}
#  \item{na.rm}{If @TRUE, non-defined unit indices are excluded, otherwise
#      not.}
# }
#
# \value{
#   Returns an @integer @vector.
# }
#
# \seealso{
#   @seemethod "getData".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getUnitIndices", "GenomeInformation", function(this, ..., na.rm=TRUE) {
  df <- getData(this, fields="chromosome", ...)
  df <- rownames(df)
  suppressWarnings({  # Suppress warnings about NAs
    df <- as.integer(df)
  })
  if (na.rm)
    df <- df[!is.na(df)]
  df
}, protected=TRUE)
