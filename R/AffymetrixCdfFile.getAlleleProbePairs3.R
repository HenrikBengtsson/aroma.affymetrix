###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod getAlleleProbePairs3
#
# @title "Gets the indices of probepairs with the same pair of SNP nucleotides"
#
# \description{
#   @get "title".
#   Note that the order of allele A and allele B is irrelevant.
#   For instance, all probepairs with nucleotides (A,T) are calibrated
#   together with all probepairs with nucleotides (T,A) reversed.
# }
#
# @synopsis
#
# \arguments{
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @list where each element is a two-column @matrix where
#   the column names are the nucleotides for the two alleles.
# }
#
# \section{Benchmarking}{
#   On an IBM Thinkpad A31 with 1.8GHz and 1GB RAM:
#   \itemize{
#    \item{Mapping10K\_Xba142}{10208 units & 432964 cells: 11 seconds.}
#    \item{Mapping50K\_Xba240}{58960 SNPs & 589,600 (PMA,PMB) probe pairs: 11 seconds.}
#   }
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAlleleProbePairs3", "AffymetrixCdfFile", function(this, units=NULL, ignoreOrder=TRUE, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Identifying the probes stratified by allele basepairs")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this)
  key <- list(method="getAlleleProbePairs", class=class(this)[1], version="2008-08-31", chipType=chipType, units=units, ignoreOrder=ignoreOrder)
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getAlleleProbePairs3", chipType=chipType, units=units, ignoreOrder=ignoreOrder)
  }
  dirs <- c("aroma.affymetrix", chipType)
  if (!force) {
    res <- loadCache(key=key, dirs=dirs)
    if (!is.null(res)) {
      verbose && cat(verbose, "Loaded from file cache")
      verbose && exit(verbose)
      return(res)
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify all possible allele pairs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading all possible allele basepairs")
  unitsWanted <- units

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify genotype units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying genotyping units (with 2 or 4 groups)")
  # Use only units that are SNPs...
  verbose && enter(verbose, "Reading unit types for all units")
  unitTypes <- getUnitTypes(this, verbose=less(verbose, 1))
  verbose && exit(verbose)

  units <- which(unitTypes == 2)
  verbose && cat(verbose, "Number of SNP units: ", length(units))
  # Not needed anymore
  unitTypes <- NULL

  # ...and with either 2 or 4 groups
  verbose && enter(verbose, "Reading number of groups per SNP unit")
  unitSizes <- nbrOfGroupsPerUnit(this, units=units)
  verbose && cat(verbose, "Detected unit sizes:")
  verbose && print(verbose, table(unitSizes))
  verbose && exit(verbose)

  units <- units[is.element(unitSizes, c(2,4))]
  # Not needed anymore
  unitSizes <- NULL
  verbose && cat(verbose, "Number of SNP units with 2 or 4 groups: ",
                                                            length(units))
  verbose && exit(verbose)

  # Clean up
  gc <- gc()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subsetting?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Operate only on a subset of probes?
  if (!is.null(unitsWanted)) {
    verbose && enter(verbose, "Subsetting")
    units <- intersect(units, unitsWanted)
    verbose && exit(verbose)
  }

  nbrOfUnits <- length(units)
  verbose && cat(verbose, "Number of SNP units to query: ", nbrOfUnits)
  if (nbrOfUnits == 0) {
    verbose && exit(verbose)
    return(NULL)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the probe nucleotides
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying probe nucleotides")

  complementaryMap <- c(A="T", C="G", G="C", T="A")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- getPathname(this)
  verbose && cat(verbose, "CDF pathname: ", pathname)

  nbrOfUnitsPerChunk <- 50e3
  nbrOfChunks <- ceiling(nbrOfUnits / nbrOfUnitsPerChunk)

  unitsTodo <- units
  count <- 1L
  allCdfData <- list()
  while (length(unitsTodo) > 0) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", count, nbrOfChunks))
    if (length(unitsTodo) <= nbrOfUnitsPerChunk) {
      uu <- 1:length(unitsTodo)
    } else {
      uu <- 1:nbrOfUnitsPerChunk
    }

    unitsChunk <- unitsTodo[uu]
    unitsTodo <- unitsTodo[-uu]
    verbose && cat(verbose, "Units: ")
    verbose && str(verbose, unitsChunk)

    verbose && enter(verbose, "Reading group names, group directions, and cell indices")
    cdfData <- .readCdf(pathname, units=unitsChunk, readIndices=TRUE, readXY=FALSE, readAtoms=FALSE, readIndexpos=FALSE, readBases=FALSE, readUnitNumber=FALSE, readUnitType=FALSE, readUnitDirection=FALSE, readGroupAtomNumbers=FALSE, stratifyBy="pm")
    verbose && exit(verbose)

    # Save memory by removing names. [55Mb -> 44Mb]
    names(cdfData) <- NULL

    verbose && enter(verbose, "Extracting fields of interest")
    cdfData <- lapply(cdfData, FUN=.subset2, "groups")
    verbose && exit(verbose)

    verbose && enter(verbose, "Pairing up cell indices")
    verbose && printf(verbose, "Progress: %d, ", length(unitsChunk))
    for (uu in seq_along(cdfData)) {
      if (uu %% 1000 == 0)
        verbose && writeRaw(verbose, length(unitsChunk)-uu, ", ")
      unit <- cdfData[[uu]]
      nucleotides <- names(unit)
      directions <- sapply(unit, .subset2, "groupdirection")
      directions <- unlist(directions, use.names=FALSE)
      swap <- which(directions == "sense")
      if (length(swap) > 0) {
        nucleotides[swap] <- complementaryMap[nucleotides[swap]]
      }
      indices <- lapply(unit, .subset2, "indices")
      nbrOfPairs <- length(indices) %/% 2
      unit <- list()
      for (gg in 1:nbrOfPairs) {
        idxs <- 2*(gg-1)+1:2
        pair <- nucleotides[idxs]
        values <- matrix(unlist(indices[idxs], use.names=FALSE), ncol=2)
        group <- gg
        if (ignoreOrder) {
          o <- order(pair)
          if (o[1] == 2) {
            pair <- pair[o]
            values <- values[,o,drop=FALSE]
            group <- -group
          }
        }
        pair <- paste(pair, collapse="")
        data <- cbind(unitsChunk[uu], group, values)
        dimnames(data) <- NULL
        data <- t(data)
        unit[[pair]] <- c(unit[[pair]], data)
      } # for (gg ...)

      cdfData[[uu]] <- unit
    } # for (uu ...)
    verbose && writeRaw(verbose, "0.\n")
    verbose && exit(verbose)

    # Flatten
##    cdfData <- allCdfData
    knownPairs <- lapply(cdfData, FUN=names)
    knownPairs <- unlist(knownPairs, use.names=FALSE)
    knownPairs <- sort(unique(knownPairs), na.last=TRUE)
    cellGroups <- vector("list", length(knownPairs))
    names(cellGroups) <- knownPairs
    for (pair in knownPairs) {
      values <- lapply(cdfData, FUN=.subset2, pair)
      values <- unlist(values, use.names=FALSE)
      cellGroups[[pair]] <- c(cellGroups[[pair]], values)
    }
    # Not needed anymore
    cdfData <- NULL

##    # Turn each group into a matrix
##    for (kk in seq_along(cellGroups)) {
##      values <- matrix(cellGroups[[kk]], nrow=4)
##      rownames(values) <- c("unit", "group", "A", "B")
##      cellGroups[[kk]] <- values
##      # Not needed anymore
##      values <- NULL
##    }

    # Append
    allCdfData <- c(allCdfData, list(cellGroups))

    gc <- gc()
    verbose && print(verbose, gc)

    count <- count + 1
    verbose && exit(verbose)
  } # while(...)
  # Not needed anymore
  # Not needed anymore
  units <- unitsTodo <- uu <- NULL
  verbose && exit(verbose)



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Merge chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Flatten
  cdfData <- allCdfData
  # Not needed anymore
  allCdfData <- NULL
  knownPairs <- lapply(cdfData, FUN=names)
  knownPairs <- unlist(knownPairs, use.names=FALSE)
  knownPairs <- sort(unique(knownPairs), na.last=TRUE)
  cellGroups <- vector("list", length(knownPairs))
  names(cellGroups) <- knownPairs
  for (pair in knownPairs) {
    values <- lapply(cdfData, FUN=.subset2, pair)
    values <- unlist(values, use.names=FALSE)
    cellGroups[[pair]] <- c(cellGroups[[pair]], values)
  }
  # Not needed anymore
  cdfData <- NULL

  # Turn each group into a matrix
  for (kk in seq_along(cellGroups)) {
    pair <- names(cellGroups)[kk]
    pair <- strsplit(pair, split="", fixed=TRUE)[[1]]
    values <- matrix(cellGroups[[kk]], nrow=4)
    rownames(values) <- c("unit", "group", pair)
    cellGroups[[kk]] <- values
    # Not needed anymore
    values <- NULL
  }

  gc <- gc()
  verbose && print(verbose, gc)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Summarize (for verbose output)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfAllelePairs <- length(cellGroups)
  verbose && printf(verbose, "Identified %d (PM_A,PM_B) allele-pair groups in %d units\n", nbrOfAllelePairs, nbrOfUnits)

  counts <- sapply(cellGroups, FUN=ncol)
  verbose && cat(verbose, "Number of probe pairs in each group:")
  verbose && print(verbose, counts)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving cell indices for all other units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying cell indices for all non-SNP units")
  unitTypes <- getUnitTypes(this, verbose=verbose);  # Takes time
  verbose && cat(verbose, "Table of identified unit types:")
  verbose && print(verbose, table(unitTypes))
  units <- which(unitTypes != 2)
  # Not needed anymore
  unitTypes <- NULL

  if (!is.null(unitsWanted)) {
    verbose && enter(verbose, "Subsetting")
    units <- intersect(units, unitsWanted)
    verbose && exit(verbose)
  }

  verbose && cat(verbose, "Identified non-SNP units:")
  verbose && str(verbose, units)
  if (length(units) == 0) {
    cells <- NULL
  } else {
    cells <- getCellIndices(this, units=units, useNames=FALSE, unlist=TRUE,
                                                          verbose=verbose)
  }
  # Not needed anymore
  units <- NULL
  verbose && cat(verbose, "Identified non-SNP cells:")
  verbose && str(verbose, cells)
  verbose && exit(verbose)


  res <- list(snps=cellGroups, nonSNPs=cells)
  # Not needed anymore
  cellGroups <- cells <- NULL


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Saving to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  saveCache(res, key=key, dirs=dirs)

  verbose && exit(verbose)

  res
}, protected=TRUE) # getAlleleProbePairs3()
