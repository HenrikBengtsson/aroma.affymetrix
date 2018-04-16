setMethodS3("getAlleleProbePairs2", "AffymetrixCdfFile", function(this, ..., verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  cdfGetGroups <- affxparser::cdfGetGroups

  # Look up base::apply(); '::' is expensive
  base_apply <- base::apply


  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)

  verbose && enter(verbose, "Identifying the probes stratified by allele basepairs")
  cdfFile <- getPathname(this)

  # Identify all possible allele pairs
  verbose && enter(verbose, "Loading all possible allele basepairs")
  # Use only units that are SNPs
  types <- getUnitTypes(this, verbose=verbose)
  units <- which(types == 2)

  # Read group names for the SNPs
  groupNames <- .readCdfGroupNames(cdfFile, units=units)
  uGroupNames <- unique(groupNames)
  verbose && exit(verbose)

  uGroupNames0 <- lapply(uGroupNames, FUN=function(x) {
    x <- matrix(x, nrow=2)
    if (ncol(x) == 2) {
      # Take the complement bases for the reverse strand
      x[,2] <- c(A="T", C="G", G="C", T="A")[x[,2]]
    }
    x

  })

  uBasepairs0 <- lapply(uGroupNames0, FUN=function(x) {
    base_apply(x, MARGIN=2, FUN=sort)
  })

  uBasepairs1 <- lapply(uBasepairs0, FUN=function(x) {
    base_apply(x, MARGIN=2, FUN=paste, collapse="")
  })

  # Get all unique allele basepairs
  uBasepairs <- sort(unique(unlist(uBasepairs1)))
  uBasepairs <- strsplit(uBasepairs, split="")

  # Create basepairs to group names map
  map <- vector("list", length(uBasepairs))
  names(map) <- uBasepairs
  bpIdx <- vector("list", length(uBasepairs0))
  for (bp in uBasepairs) {
    for (kk in 1:length(bpIdx)) {
      set <- uBasepairs0[[kk]]
      bpIdx[[kk]] <- kk + which(bp == set)/10
    }
    map[[bp]] <- unlist(bpIdx)
  }

  # Read all of the CDF file
  verbose && enter(verbose, "Loading cell indices for all probepairs")
  cdfAll <- .readCdfCellIndices(cdfFile, units=units, stratifyBy="pm")
  # Not needed anymore
  units <- NULL
  verbose && exit(verbose)

  verbose && enter(verbose, "Stratifying by unique allele basepairs")
  probes <- vector("list", length(map))
  for (kk in 1:length(map)) {
    basepair <- names(map)[kk]
    verbose && enter(verbose, "Allele basepair ", basepair)

    bpIdx <- map[[kk]]
    gIdx <- as.integer(bpIdx)
    sIdx <- round(10*(bpIdx - gIdx))

    verbose && cat(verbose, "Located in ", length(unique(gIdx)), " group(s).")

    idx <- lapply(groupNames, FUN=identical, basepair)
    idx <- which(unlist(idx, use.names=FALSE))
    cdf <- cdfAll[idx]

    cdf0 <- vector("list", length=4)
    for (gg in 1:4) {
      cells <- .applyCdfGroups(cdf, cdfGetGroups, gg)
      cells <- unlist(cells, use.names=FALSE)
      cdf0[[gg]] <- cells
    }
    probes[[kk]] <- cdf0
    names(probes)[kk] <- basepair

    verbose && exit(verbose)
  }
  verbose && exit(verbose)

  probes
}, private=TRUE) # getAlleleProbePairs2()
