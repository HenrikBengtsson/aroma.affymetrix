# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: Platform specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("importFromGenomeInformation", "AromaUgpFile", function(this, gi, ..., verbose=FALSE) {
  # Argument 'gi':
  gi <- Arguments$getInstanceOf(gi, "GenomeInformation")

  # AD HOC patch, since units==NULL does not work./HB 2007-03-03
  units <- seq_len(nbrOfUnits(gi))
  data <- getData(gi, units=units, fields=c("chromosome", "physicalPosition"))

  chr <- data[,"chromosome"]
  if (is.character(chr)) {
    chr[chr == "X"] <- 23
    chr[chr == "Y"] <- 24
    chr[chr %in% c("MT", "Z")] <- 25
    suppressWarnings({
      chr <- as.integer(chr)
    })
  }
  
  pos <- data[,"physicalPosition"]
  suppressWarnings({
    pos <- as.integer(pos)
  })

  this[,1] <- chr
  this[,2] <- pos

  # A best guess of what was imported
  units <- units[!(is.na(chr) & is.na(pos))]

  invisible(units)
})



setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUgpFile", function(this, csv, shift="auto", onReplicates=c("median", "mean", "overwrite"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csv':
  csv <- Arguments$getInstanceOf(csv, "AffymetrixNetAffxCsvFile")

  # Argument 'shift':
  if (!identical(shift, "auto"))
    shift <- Arguments$getInteger(shift)

  # Argument 'onReplicates':
  onReplicates <- match.arg(onReplicates)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Importing (unit name, chromosome, position) data from ", class(csv)[1])

  # Query CDF
  unf <- getUnitNamesFile(this)
  unfUnitNames <- getUnitNames(unf)

  # Read data
  data <- readDataUnitChromosomePosition(csv, ..., verbose=less(verbose))
  importNames <- attr(data, "importNames")

  # Map to unit names
  unfUnits <- match(data[,1], unfUnitNames)

  # Exclude units that are not in the annotation unit names file
  keep <- which(!is.na(unfUnits))
  unfUnits <- unfUnits[keep]
  if (length(unfUnits) == 0) {
    warning("None of the imported unit names match the ones in the annotation unit names file ('", getPathname(unf), "'). Is the correct file ('", getPathname(csv), "'), being imported?")
  }
  data <- data[keep,2:3,drop=FALSE]
  importNames <- importNames[2:3]

  # Garbage collect
  # Not needed anymore
  keep <- NULL
  gc <- gc()


  # Shift positions?
  if (identical(shift, "auto")) {
    shift <- 0
    if (any(regexpr("[sS]tart", importNames) != -1))
      shift <- 13
  }
  if (shift != 0) {
    verbose && printf(verbose, "Shifting positions %d steps.\n", shift)
    data[,2] <- data[,2] + as.integer(shift)
  }
 
  # Replicated positions per unit?
  dups <- which(duplicated(unfUnits))
  if (length(dups) > 0) {
    verbose && enter(verbose, "Detected units with replicated positions")
    dupUnits <- unique(unfUnits[dups])
    nDupUnits <- length(dupUnits)
    verbose && cat(verbose, "Number of units with replicated positions: ", 
                                                               nDupUnits)

    if (onReplicates %in% c("median", "mean")) {
      verbose && enter(verbose, "Calculate average positions for those (assuming they are on the same chromosome)")
      if (onReplicates == "mean") {
        avgFcn <- median
      } else {
        avgFcn <- mean
      }
      for (kk in seq_along(dupUnits)) {
        if (kk %% 500 == 0)
         verbose && printf(verbose, "%d, ", kk)
        dupUnit <- dupUnits[kk]
        # Identify position
        units <- which(unfUnits == dupUnit)
        # Average position
        avgPos <- median(data[units,2], na.rm=TRUE)
        avgPos <- round(avgPos)
        # Update (can we update just units[1]?)
        data[units,2] <- avgPos
      }
      verbose && cat(verbose, kk)
      verbose && exit(verbose)
      verbose && enter(verbose, "Remove the extraneous cases")
      data <- data[-dups,,drop=FALSE]
      unfUnits <- unfUnits[-dups]
      verbose && exit(verbose)
      # Not needed anymore
      dupUnits <- units <- NULL
      verbose && str(verbose, unfUnits)
      warning("The positions for ", nDupUnits, " units were calculated as the average of replicated positions, since that was what was available on file.")
    } else {
      verbose && cat(verbose, "Ignored replicated positions. The last one written will be the one available.")
    }
    verbose && exit(verbose)
  }
  # Not needed anymore
  dups <- NULL

  # Garbage collect
  gc <- gc()
  verbose && print(verbose, gc, level=-10)

  this[unfUnits,1] <- data[,1]
  this[unfUnits,2] <- data[,2]

  # Not needed anymore
  data <- NULL
  gc <- gc()
  verbose && print(verbose, gc, level=-10)

  verbose && exit(verbose)

  invisible(unfUnits)
})



setMethodS3("importFromAffymetrixTabularFile", "AromaUgpFile", function(this, src, colClasses=c("*"="NULL", "^probeSetID$"="character", "^chromosome$"="character", "^(physicalPosition|position)$"="character"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'src':
  src <- Arguments$getInstanceOf(src, "AffymetrixTabularFile")

  units <- importFromGenericTabularFile(this, src=src, 
            colClasses=colClasses, camelCaseNames=TRUE, ...)

  invisible(units)
}, protected=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: Platform specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
