setMethodS3("getSetsOfProbes", "AllelicCrosstalkCalibration", function(this, ..., version=c(0,1,3,4), fakeSymmetry=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'version':
  version <- Arguments$getIntegers(version)
  version <- version[1]
  if (version == 0) {
    if (this$.pairBy == "CDF") {
      version <- 1
    } else if (this$.pairBy == "sequence") {
      version <- 4
    } else {
      version <- 4
    }
  }

  # Argument 'fakeSymmetry':
  fakeSymmetry <- Arguments$getLogical(fakeSymmetry)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  setsOfProbes <- this$.setsOfProbes

  if (force || is.null(setsOfProbes)) {
    verbose && enter(verbose, "Identifying sets of pairs of cell indices")
    dataSet <- getInputDataSet(this)
    cdf <- getCdf(dataSet)
    chipType <- getChipType(cdf)
    verbose && cat(verbose, "Chip type: ", chipType)

    params <- getParameters(this)
    mergeShifts <- params$mergeShifts
    verbose && cat(verbose, "Merge shifts: ", mergeShifts)
    B <- params$B
    verbose && cat(verbose, "Number of nucleotides: ", B)

    # Assert that needed annotation files exist
    needAcs <- FALSE
    if (version == 4) {
      needAcs <- TRUE
    } else if (B > 1 || !mergeShifts) {
      needAcs <- TRUE
    }

    if (needAcs) {
      verbose && enter(verbose, "Locating AromaCellSequenceFile")
      chipType <- getChipType(cdf, fullname=FALSE)
      nbrOfCells <- nbrOfCells(cdf)
      verbose && cat(verbose, "Chip type: ", chipType)
      verbose && cat(verbose, "Number of cells: ", nbrOfCells)
      acs <- AromaCellSequenceFile$byChipType(chipType, nbrOfCells=nbrOfCells)
      verbose && print(verbose, acs)
      verbose && exit(verbose)
    } else {
      acs <- NULL
    }


    if (version == 1) {
      setsOfProbes <- getAlleleProbePairs(cdf, verbose=verbose)
      # Coerce to new structure. /HB 2008-09-02
      nonSNPs <- setsOfProbes$nonSNPs
      setsOfProbes$nonSNPs <- NULL
      setsOfProbes <- lapply(setsOfProbes, FUN=t)
      pairs <- names(setsOfProbes)
      pairs <- strsplit(pairs, split="", fixed=TRUE)
      pairs <- sapply(pairs, FUN=paste, collapse="/")
      names(setsOfProbes) <- pairs
      setsOfProbes <- list(snps=setsOfProbes, nonSNPs=nonSNPs)
      # Not needed anymore
      nonSNPs <- pairs <- NULL
    } else if (version == 3) {
      setsOfProbes <- getAlleleProbePairs3(cdf, ignoreOrder=TRUE, verbose=verbose)
      snps <- setsOfProbes$snps
      for (kk in seq_along(snps)) {
        cells <- snps[[kk]][3:4,,drop=FALSE]
        o <- order(cells[1,])
        cells <- cells[,o,drop=FALSE]
        snps[[kk]] <- cells
        # Not needed anymore
        cells <- o <- NULL
      }
      pairs <- names(snps)
      pairs <- strsplit(pairs, split="", fixed=TRUE)
      pairs <- sapply(pairs, FUN=paste, collapse="/")
      names(snps) <- pairs
      setsOfProbes$snps <- snps
      # Not needed anymore
      snps <- pairs <- NULL
    } else if (version == 4) {
      verbose && enter(verbose, "Identifying cell indices for all non-SNP units")
      unitTypes <- getUnitTypes(cdf, verbose=verbose)
      units <- which(unitTypes != 2)
      # Not needed anymore
      unitTypes <- NULL
      verbose && enter(verbose, "Non-SNP units:")
      verbose && str(verbose, units)
      if (length(units) > 0) {
        nonSNPs <- getCellIndices(cdf, units=units,
                       useNames=FALSE, unlist=TRUE, verbose=verbose)
      } else {
        nonSNPs <- NULL
      }
      # Not needed anymore
      units <- NULL
      verbose && enter(verbose, "Non-SNP cells:")
      verbose && str(verbose, nonSNPs)
      verbose && exit(verbose)

      cells <- getAlleleCellPairs(cdf, verbose=verbose)
      S <- as.integer(B %/% 2 + 1)
      shifts <- -(S-1):(S-1)
      verbose && cat(verbose, "Probe shifts:")
      verbose && print(verbose, shifts)
      snps <- groupBySnpNucleotides(acs, cells=cells, shifts=shifts,
                                                     verbose=verbose)
      # Not needed anymore
      cells <- shifts <- NULL
      # Clean out empty sets
      for (kk in seq_along(snps)) {
        cells <- snps[[kk]]
        if (length(cells) == 0)
          snps[[kk]] <- NULL
      }
      setsOfProbes <- list(snps=snps, nonSNPs=nonSNPs)
      # Not needed anymore
      snps <- nonSNPs <- NULL
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Group by number of nucleotides (B) around SNP position?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (B == 0) {
      verbose && enter(verbose, "Merging SNP groups")
      snpsT <- list(all=NULL)
      snps <- setsOfProbes$snps
      for (kk in seq_along(snps)) {
        snpsT$all <- cbind(snpsT$all, snps[[kk]])
      }
      # Not needed anymore
      snps <- NULL
      setsOfProbes$snps <- snpsT
      # Not needed anymore
      snpsT <- NULL
      verbose && exit(verbose)
    } else if (B == 1) {
      # Nothing to do, default
    } else {
      if (version != 4)
        throw("Not supported.")
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Merge shifts?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (mergeShifts) {
      # Nothing to do.
    } else {
      # Split by shift
      verbose && enter(verbose, "Grouping probe pairs by their shift relative to SNP")

      # For each set of cells...
      for (gg in seq_along(setsOfProbes)) {
        cells <- setsOfProbes[[gg]]
        groupTag <- names(setsOfProbes)[gg]
        verbose && enter(verbose, sprintf("Set #%d ('%s') of %d", gg, groupTag, length(setsOfProbes)))

        dim <- dim(cells)
        if (is.null(dim) || dim[2] != 2) {
          # Not a set of SNP cell indexes.
          next
        }

        # Identify the location of the SNP position for each probe pair
        snpPositions <- rep(NA_integer_, times=dim[1])
        possibleShifts <- as.integer(seq(from=-4, to=+4));  # Hardwired!
        possiblePositions <- 13 + possibleShifts
        for (pos in possiblePositions) {
          seqsA <- readSequenceMatrix(acs, cells=cells[,1], positions=pos)
          seqsB <- readSequenceMatrix(acs, cells=cells[,2], positions=pos)
          snpPositions[(seqsA != seqsB)] <- pos
          # Not needed anymore
          seqsA <- seqsB <- NULL
        }
        shifts <- snpPositions - as.integer(13)

        possibleShifts <- sort(unique(shifts), na.last=TRUE)

        # For each shift...
        subgroups <- vector("list", length(possibleShifts))
        shiftTags <- sprintf("shift=%+d", possibleShifts)
        names(subgroups) <- paste(groupTag, shiftTags, sep=",")

        for (ss in seq_along(possibleShifts)) {
          shift <- possibleShifts[ss]
          if (is.na(shift)) {
            idxs <- which(is.na(shifts))
          } else {
            idxs <- which(shifts == shift)
          }
          cellsSS <- cells[idxs,,drop=FALSE]
          subgroups[[ss]] <- cellsSS
          # Not needed anymore
          idxs <- cellsSS <- NULL
        } # for (shift ...)
        # Not needed anymore
        cells <- NULL

        setsOfProbes[[gg]] <- subgroups
        # Not needed anymore
        subgroups <- NULL

        verbose && exit(verbose)
      } # for (gg ...)
      verbose && exit(verbose)
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fake symmetry by flipping every 2nd (A,B) to (B,A)?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (fakeSymmetry) {
      verbose && enter(verbose, "Flipping every 2nd (A,B) to (B,A)")
      snps <- setsOfProbes$snps
      for (kk in seq_along(snps)) {
        idxs <- snps[[kk]]
        cc <- seq(from=1, to=ncol(idxs), by=2)
        idxs[,cc] <- idxs[2:1,cc, drop=FALSE]
        snps[[kk]] <- idxs
        # Not needed anymore
        cc <- idxs <- NULL
      } # for (gg ...)
      setsOfProbes$snps <- snps
      # Not needed anymore
      snps <- NULL
      verbose && exit(verbose)
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Transpose
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # For each set of cells...
    snps <- setsOfProbes$snps
    for (kk in seq_along(snps)) {
      snps[[kk]] <- t(snps[[kk]])
    } # for (gg ...)
    setsOfProbes$snps <- snps
    # Not needed anymore
    snps <- NULL

    # Tag the result with a version number
    attr(setsOfProbes, "version") <- version

    this$.setsOfProbes <- setsOfProbes

    verbose && exit(verbose)
  }

  setsOfProbes
}, protected=TRUE) # getSetsOfProbes()
