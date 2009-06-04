setMethodS3("getProbeSequenceData", "AffymetrixCdfFile", function(this, paths=NULL, rows=NULL, safe=TRUE, force=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'safe':
  safe <- Arguments$getLogical(safe);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving probe-sequence data");
  chipTypeFull <- getChipType(this, fullname=TRUE);
  verbose && cat(verbose, "Chip type (full): ", chipTypeFull);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate find probe sequence file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Locating probe-tab file");
  # The probe-sequence does not depend on the CDF but only the chip type,
  # which is why we ignore any tags for the CDF.
  chipType <- getChipType(this, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);

  ptf <- AffymetrixProbeTabFile$byChipType(chipType=chipType, 
                                                 verbose=less(verbose, 100));
  verbose && print(verbose, ptf);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate content against CDF?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (safe) {
    verbose && enter(verbose, "Validating probe-tab file against CDF");

    # Reading the first unit name
    data <- readDataFrame(ptf, colClassPatterns=c("^unitName$"="character"), 
                                          rows=1, verbose=less(verbose, 50));
    unitName <- data$unitName;
    verbose && str(verbose, "Unit name: ", unitName);
    unit <- indexOf(this, names=unitName);
    verbose && cat(verbose, "Unit index: ", unit);
    if (is.na(unit)) {
      throw("Failed to identify CDF unit with unit name '", unitName, "': ",
                                                         getPathname(ptf));
    }

    # Reading the (x,y) & sequence data
    data <- readSequenceDataFrame(ptf, rows=1, verbose=less(verbose, 50));
    verbose && print(verbose, data);

    xySeq <- c(data$probeXPos, data$probeYPos);
    verbose && cat(verbose, "(x,y):");
    verbose && print(verbose, xySeq);

    unitInfo <- readUnits(this, units=unit);
    x <- applyCdfGroups(unitInfo, cdfGetFields, "x");
    x <- unlist(x, use.names=FALSE);
    y <- applyCdfGroups(unitInfo, cdfGetFields, "y");
    y <- unlist(y, use.names=FALSE);

    # Now, find that (x,y) coordinate in the CDF file
    idx <- whichVector(xySeq[1] == x & xySeq[2] == y);
    # Sanity check
    if (length(idx) != 1) {
      throw("The (x,y) coordinate (", paste(xySeq, collapse=","), ") of the probe-tab file could not be found in CDF unit #", unit, " ('", unitName, "'): ", getPathname(ptf));
    }

    verbose && exit(verbose);
  } # if (safe)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading probe sequence data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading (x,y,sequence) data");
  data <- readSequenceDataFrame(ptf, rows=rows, verbose=less(verbose, 20));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating (x,y)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Validating (x,y) against CDF dimension");
  # The dimension of the chip type
  dimension <- getDimension(this);
  verbose && cat(verbose, "CDF dimension:");
  verbose && print(verbose, dimension);

  # Sanity check
  xRange <- sprintf("[0,%d]", dimension[1]);
  yRange <- sprintf("[0,%d]", dimension[2]);

  if (any(data$probeXPos < 0 | data$probeXPos > dimension[1]-1)) {
    throw("Detected probe x position out of range ", xRange, ": ",
                                                getPathname(ptf));
  }
  if (any(data$probeYPos < 0 | data$probeYPos > dimension[2]-1)) {
    throw("Detected probe y position out of range ", yRange, ": ",
                                                 getPathname(ptf));
  }
  verbose && exit(verbose);
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reformating data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Renaming columns");
  names <- colnames(data);
  names[names == "probeXPos"] <- "x";
  names[names == "probeYPos"] <- "y";
  names[names == "probeSequence"] <- "sequence";
  colnames(data) <- names;
  verbose && exit(verbose);

  verbose && enter(verbose, "Calculating and appending cell indices");
  dimension <- getDimension(this);
  cells <- data$y * dimension[1] + data$x + as.integer(1);
  data <- cbind(cell=cells, data);
  verbose && str(verbose, data);
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  data;
}, private=TRUE)



###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod computeAffinities
#
# @title "Calculates probe affinities from sequence"
#
# \description{
#  @get "title".
#
# Adapted from @see "gcrma::compute.affinities" in the \pkg{gcrma} package.
# Attempts to find the tab-separated probe sequence file associated with
# a particular CDF, and matches sequence to probe index in order to assign
# an affinity to each probe.
# }
#
# @synopsis
#
# \arguments{
#   \item{paths}{A @character variable containing location(s) to look for
#    the probe sequence file.  Multiple locations should be separated by
#    semi-colons.}
#   \item{force}{If @FALSE, cached results is returned, if available.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @numeric @vector of (log2) probe affinities, of length equal
#  to the total number of features on the array.
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#*/###########################################################################
setMethodS3("computeAffinities", "AffymetrixCdfFile", function(this, paths=NULL, safe=TRUE, force=FALSE, verbose=FALSE, ...) {
  isPackageInstalled("gcrma") || throw("Package not installed: gcrma");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # getSeqMatrix() corresponds to calling the following native function
  #   .Call("gcrma_getSeq", seq, PACKAGE="gcrma");
  # However, we do want to avoid to load 'gcrma', because it loads Biostrings
  # which loads IRanges, which cause problems, e.g. trim(). /HB 2009-05-09
  getSeqMatrix <- function(seq, ...) {
    seq <- charToRaw(seq);
    map <- charToRaw("ACGT");
    res <- matrix(0L, nrow=4, ncol=length(seq));
    for (kk in 1:4) {
      res[kk,(seq == map[kk])] <- 1L;
    }
    res;
  } # getSeqMatrix()

  getGcrmaSplineCoefs <- function(...) {
    env <- new.env();
    data("affinity.spline.coefs", package="gcrma", envir=env);
    env$affinity.spline.coefs;
  } # getGcrmaSplineCoefs()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'safe':
  safe <- Arguments$getLogical(safe);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  chipTypeFull <- getChipType(this, fullname=TRUE);


  verbose && enter(verbose, "Computing GCRMA probe affinities for ", nbrOfUnits(this), " units");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Looking for MMs (and PMs) in the CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identify PMs and MMs among the CDF cell indices");
  isPm <- isPm(this);
  verbose && str(verbose, isPm);
  verbose && summary(verbose, isPm);
  nbrOfPMs <- sum(isPm);
  verbose && cat(verbose, "MMs are defined as non-PMs");
  isMm <- !isPm;
  nbrOfMMs <- sum(isMm);
  verbose && cat(verbose, "Number of PMs: ", nbrOfPMs);
  verbose && cat(verbose, "Number of MMs: ", nbrOfMMs);
  verbose && exit(verbose);

  # Sanity check
  if (nbrOfMMs == 0) {
#    throw("Cannot calculate gcRMA probe affinities. The CDF contains no MMs: ", chipTypeFull);
  }

  # Checking cache
  key <- list(method="computeAffinities", class=class(this)[1], 
              chipTypeFull=chipTypeFull, version="2009-05-09");
  dirs <- c("aroma.affymetrix", chipTypeFull);
  if (!force) {
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(res);
    }
  }


  verbose && enter(verbose, "Reading probe-sequence data");
  sequenceInfo <- getProbeSequenceData(this, safe=safe, verbose=verbose);
  verbose && str(verbose, sequenceInfo);
  verbose && exit(verbose);

  nbrOfSequences <- nrow(sequenceInfo);

  verbose && enter(verbose, "Get CDF cell indices for all units");
  cells <- getCellIndices(this, useNames=FALSE, unlist=TRUE);
  verbose && exit(verbose);

  # probe sequence tab-separated files generally only contain the PM
  # sequence (since we know that MM sequence always has base 13 complemented).
  # This is true for expression arrays, but possibly not for genotyping arrays.
  # We will calculate affinity for PM and MM anyway, and then
  # try to match MM indices from the CDF file with their appropriate PMs (and
  # hence assign the correct affinity).
   
  # assume that MM has same x-coordinate as PM, but y-coordinate larger by 1
  indexPm <- sequenceInfo$cell;
  dimension <- getDimension(this);
  indexMmPutative <- indexPm + dimension[1];

  # match up putative MM with actual list of MMs from CDF
  matches <- match(indexMmPutative, cells[isMm]);
  indexMm <- cells[isMm][matches];
  notNA <- which(!is.na(indexMm));
  indexMm <- indexMm[notNA];
  
  # for PM+MM arrays, the number of MMs and the number of PMs in the
  # CDF is equal
  isPMMMChip <- (length(indexMm) == length(indexPm));
  
  rm(indexMmPutative, matches); # Not needed anymore


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate basis matrix based on priori probe affinitites (of gcrma)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # now calculate affinities - code reused from compute.affinities() in gcrma
  verbose && enter(verbose, "Calculating probe affinities");

  # To please R CMD check on R v2.6.0
  affinity.spline.coefs <- getGcrmaSplineCoefs();

  affinity.basis.matrix <- splines::ns(1:25, df=length(affinity.spline.coefs)/3);

  A13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[1:5]);
  T13 <- 0;
  C13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[6:10]);
  G13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[11:15]);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Predicting probe affinities based on the probe sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfCells <- dimension[1]*dimension[2];
  rm(dimension); # Not needed anymore

  # Allocate empty vector of affinities
  naValue <- as.double(NA);
  affinities <- rep(naValue, nbrOfCells);

  if(isPMMMChip) {
    # ---------------------------------------------------------
    # keep this code below untouched for the PM+MM expression arrays
    # ---------------------------------------------------------
    apm <- vector("double", nbrOfSequences);
    amm <- vector("double", nbrOfSequences);

    T13A13 <- T13 - A13;
    C13G13 <- C13 - G13;

    sequences <- sequenceInfo$sequence;
    rm(sequenceInfo);

    if (verbose) {
      cat(verbose, "Progress (counting to 100): ");
      pb <- ProgressBar(stepLength=100/(nbrOfSequences/1000));
      reset(pb);
    }
  
    for (ii in seq(along=apm)) {
      if (verbose && ii %% 1000 == 0)
        increase(pb);

      # Get a 4x25 matrix with rows A, C, G, and T.
      charMtrx <- getSeqMatrix(sequences[ii]);

      # Calculate the PM affinity
      A <- cbind(charMtrx[1,,drop=TRUE] %*% affinity.basis.matrix, 
                 charMtrx[2,,drop=TRUE] %*% affinity.basis.matrix, 
                 charMtrx[3,,drop=TRUE] %*% affinity.basis.matrix);
      apm[ii] <- A %*% affinity.spline.coefs;

      # Calculate the MM affinity
      base <- whichVector(charMtrx[,13,drop=TRUE] == 1);
      if (base == 1) {
        amm[ii] <- apm[ii] + T13A13;  # + T13 - A13
      } else if (base == 4) {
        amm[ii] <- apm[ii] - T13A13;  # + A13 - T13
      } else if (base == 3) {
        amm[ii] <- apm[ii] + C13G13;  # + C13 - G13
      } else if (base == 2) {
        amm[ii] <- apm[ii] - C13G13;  # + G13 - C13
      } else {
        throw("Unknown nucleotide index: ", base);
      }
    } # for (ii in ...)
    rm(charMtrx, A); # Not needed anymore
    reset(pb);
    verbose && cat(verbose, "");

    # create a vector to hold affinities and assign values to the 
    # appropriate location in the vector
    affinities[indexPm] <- apm;
    affinities[indexMm] <- amm[notNA];
    rm(indexPm, indexMm, apm, amm, notNA); # Not needed anymore
    # ---------------------------------------------------------
  } else {
    # ---------------------------------------------------------
    # new code to compute affinities from the MM probes
    # (antigenomic or whatever) on PM-only arrays  -- MR 2009-03-28
    # ---------------------------------------------------------
    naValue <- as.double(NA);
    apm <- rep(naValue, times=nbrOfSequences);

    indexAll <- sequenceInfo$cell;
    sequences <- sequenceInfo$sequence;
    rm(sequenceInfo);

    # Process only sequences of length 25
    n <- nchar(sequences);
    idxs <- whichVector(n == 25);
    nbrOfNon25mers <- nbrOfSequences - length(idxs);
    if (nbrOfNon25mers > 0) {
      apm[-idxs] <- naValue;
      warning("Detected ", nbrOfNon25mers, " sequence that are not of length 25 nucleotides. For those probes, the affinities are defined to be NA.");
    }

    if (verbose) {
      cat(verbose, "Progress (counting to 100): ");
      pb <- ProgressBar(stepLength=100/(length(idxs)/1000));
      reset(pb);
    }
  
    for (ii in seq(along=idxs)) {
      if (verbose && ii %% 1000 == 0)
        increase(pb);
      idx <- idxs[ii];
      # Get a 4x25 matrix with rows A, C, G, and T.
      charMtrx <- getSeqMatrix(sequences[idx]);

      # Calculate the PM affinity
      A <- cbind(charMtrx[1,,drop=TRUE] %*% affinity.basis.matrix, 
                 charMtrx[2,,drop=TRUE] %*% affinity.basis.matrix, 
                 charMtrx[3,,drop=TRUE] %*% affinity.basis.matrix);
      apm[idx] <- A %*% affinity.spline.coefs;
    } # for (ii in ...)
    reset(pb);

    affinities[indexAll] <- apm;

    rm(indexAll, apm); # Not needed anymore
  } # if(isPMMMChip)
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Saving to cache
  comment <- paste(unlist(key, use.names=FALSE), collapse=";");
  saveCache(key=key, affinities, comment=comment, dirs=dirs);

  verbose && exit(verbose);
  
  affinities;
}, private=TRUE)

############################################################################
# HISTORY:
# 2009-05-09 [HB]
# o CLEAN UP: computeAffinities() of AffymetrixCdfFile no longer need to
#   load 'gcrma'.  However, it still needs to load a small set of 
#   prior estimates from the 'gcrma' package.
#   'matchprobes' package, but never used it.
# o CLEAN UP: computeAffinities() of AffymetrixCdfFile did load the
#   'matchprobes' package, but never used it.
# o Added private getProbeSequenceData() to AffymetrixCdfFile. This method
#   will later retrieve probe sequences from the binary ACS file instead
#   of the probe-tab files.
# o Updated to make full use of the AffymetrixProbeTabFile class, which
#   translateColumnNames() will take care of all variations of names for
#   the column containing unit names.
# 2009-03-28 [MR]
# o Made several modifications to allow computing of affinities for 
#   Gene 1.0 ST arrays.  For example:
#   -- left the PM+MM array code mostly untouched
#   -- fixed some assumptions about the columns of the probe_tab file
#   -- added a different stream for PM-only (with NCs) 
# 2007-07-30
# o UPDATE: Now computeAffinities() for AffymetrixCdfFile gives an error
#   if there are no MMs in the CDF.
# o BUG FIX: computeAffinities() for AffymetrixCdfFile searched for the
#   probe-tab file using the chip type given by the fullname of the CDF
#   and not the basic name.
# 2007-09-06
# o Made computeAffinities() more memory efficient since it is using the
#   new getCellIndices(cdf, useNames=FALSE, unlist=TRUE).
# 2007-06-07
# o BUG FIX: If an Affymetrix probe tab file is not found for the chip type,
#   computeAffinitities() of AffymetrixCdfFile would throw "Error in 
#   paste(..., sep = sep, collapse = collapse): object "pattern" not found"
#   instead of an intended and more information error.
# 2007-02-06
# o Now computeAffinities() is locating the Affymetrix' probe-sequence file
#   using AffymetrixProbeTabFile$findByChipType().
# 2006-10-09
# o Added caching to file.
# o Now default probe affinity is NA (for non estimated probes).
# o Added a progress bar to the calculations.
# o Made internal reader smart so it guess the X, Y, and sequence columns
#   by reading the data and comparing to CDF information.
# o Added verbose statements.
# o Made sure connection is closed regardless how the method exits.
# 2006-10-04
# o Debugged, tested, docs/comments added.
# 2006-10-01
# o Created.
############################################################################
