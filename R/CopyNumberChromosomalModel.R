###########################################################################/**
# @RdocClass CopyNumberChromosomalModel
#
# @title "The CopyNumberChromosomalModel class"
#
# \description{
#  @classhierarchy
#
#  This \emph{abstract} class represents a copy-number model.
# }
# 
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "CopyNumberDataSetTuple".}
#   \item{refTuple}{An optional @see "CopyNumberDataFile", 
#      or @see "CopyNumberDataSet" or @see "CopyNumberDataSetTuple" 
#      for pairwise comparisons.}
#   \item{tags}{A @character @vector of tags.}
#   \item{genome}{A @character string specifying what genome is process.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Requirements}{
#   This class requires genome information annotation files for 
#   every chip type.
# }
#
# @author
#*/###########################################################################
setConstructorS3("CopyNumberChromosomalModel", function(cesTuple=NULL, refTuple=NULL, tags="*", genome="Human", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cesTuple':
  if (!is.null(cesTuple)) {
    # Coerce to CopyNumberDataSetTuple, if needed
    cesTuple <- as.CopyNumberDataSetTuple(cesTuple);

    # Currently only total copy-number estimates are accepted
    if (hasAlleleBFractions(cesTuple)) {
      throw("Unsupported copy-number data. Currently only total copy-number estimates are supported: ", getFullName(cesTuple));
    }
  }

  # Argument 'refTuple':
  if (!is.null(refTuple)) {
    # Is the reference a single file?
    if (inherits(refTuple, "AromaMicroarrayDataFile")) {
      refTuple <- list(refTuple);
    }

    # Coerce list
    if (is.list(refTuple)) {
      refList <- refTuple;

      cesList <- getListOfSets(cesTuple);

      # Assert the number of chip types
      if (length(refList) != length(cesList)) {
        throw("The number of chip types in the references (argument 'refTuple') does not match the number of chip types in the sample (argument 'cesTuple'): ", length(refList), " != ", length(cesList));
      }

      # Coerce single reference files into reference sets
      for (kk in seq(along=refList)) {
        ref <- refList[[kk]];

        # If a data file...
        if (inherits(ref, "AromaMicroarrayDataFile")) {
          chipType <- getChipType(ref, fullname=FALSE);
          ces <- cesList[[chipType]];

          # Sanity check
          if (is.null(ces)) {
            throw("The reference (argument 'refTuple') uses a chip type not used in the sample (argument 'cesTuple'): ", chipType, " not in (", paste(names(ces), collapse=", "), ")");
          }

          # Create a data set holding a sequence of one replicated reference file
          refFiles <- rep(list(ref), nbrOfArrays(ces));
          refSet <- newInstance(ces, refFiles);
          rm(refFiles);
         
          refList[[kk]] <- refSet;
          rm(refSet);
        }
        rm(ref);
      } # for (kk ...)

      refTuple <- refList;
      rm(refList, cesList);
    } # if (is.list(refTuple))

    # Coerce to CopyNumberDataSetTuple, if needed
    refTuple <- as.CopyNumberDataSetTuple(refTuple);

    # Assert the same number of chip types in test and reference set
    if (!identical(getChipTypes(refTuple), getChipTypes(cesTuple))) {
      throw("The reference tuple has a different set of chip types compared with the test tuple");
    }

    # Assert that the reference data is of the same format as the sample data
    if (hasAlleleBFractions(refTuple) != hasAlleleBFractions(cesTuple)) {
      throw("The reference data (argument 'refTuple') is not compatible with the sample data (argument 'cesTuple'). One provides total copy numbers and the other does not.");
    }
    if (hasStrandiness(refTuple) != hasStrandiness(cesTuple)) {
      throw("The reference data (argument 'refTuple') is not compatible with the sample data (argument 'cesTuple'). One provides strand-specific data and the other does not.");
    }


    # Validate consistency between the data sets and the reference files
    cesList <- getListOfSets(cesTuple);
    refList <- getListOfSets(refTuple);
    for (kk in seq(along=cesList)) {
      ces <- cesList[[kk]];
      ref <- refList[[kk]];

      # Assert that the reference and the sample sets are of the same size
      if (nbrOfArrays(ref) != nbrOfArrays(ces)) {
        throw("The number of reference files does not match the number of sample files: ", nbrOfArrays(ref), " != ", nbrOfArrays(ces));
      }

      # Assert that the reference files are compatible with the test files
      for (jj in seq(length=nbrOfArrays(ces))) {
        cf <- getFile(ces, jj);
        rf <- getFile(ref, jj);
        if (!inherits(rf, class(cf)[1])) {
          throw(class(ref)[1], " #", kk, " of argument 'refTuple' contains file #", jj, ", that is not of the same class as the paired test file: ", class(rf)[1], " !inherits from ", class(cf)[1]);
        }
      } # for (jj ...)
    } # for (kk in ...)
  }

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }


  this <- extend(ChromosomalModel(), "CopyNumberChromosomalModel",
    .cesTuple = cesTuple,
    .refTuple = refTuple,
    .paired = (length(refTuple) > 0),
    .chromosomes = NULL,
    .tags = tags,
    .genome = genome
  );

  # Validate?
  if (!is.null(this$.cesTuple)) {
    # Validate genome
    pathname <- getGenomeFile(this);
  }

  this;
})


setMethodS3("as.character", "CopyNumberChromosomalModel", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, paste("Chip type (virtual):", getChipType(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  chipTypes <- getChipTypes(this);
  nbrOfChipTypes <- length(chipTypes);
  s <- c(s, sprintf("Number of chip types: %d", nbrOfChipTypes));
  s <- c(s, "Sample & reference file pairs:");
  cesList <- getListOfSets(getSetTuple(this));
  refList <- getRefSetTuple(this);
  if (!is.null(refList))
    refList <- getListOfSets(refList);
  for (kk in seq(along=cesList)) {
    s <- c(s, sprintf("Chip type #%d of %d ('%s'):", kk, nbrOfChipTypes, chipTypes[kk]));
    s <- c(s, "Sample data set:");
    ces <- cesList[[kk]];
    ref <- refList[[kk]];
    s <- c(s, as.character(ces));
    s <- c(s, "Reference data set/file:");
    if (is.null(ref)) {
      s <- c(s, "<average across arrays>");
    } else {
      s <- c(s, as.character(ref));
    }
  }
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, protected=TRUE)



setMethodS3("clearCache", "CopyNumberChromosomalModel", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".extractRawCopyNumbersCache")) {
    this[[ff]] <- NULL;
  }

  if (!is.null(this$.cesTuple))
    clearCache(this$.cesTuple);
  if (!is.null(this$.refTuple))
    clearCache(this$.refTuple);

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
})




setMethodS3("getRefSetTuple", "CopyNumberChromosomalModel", function(this, ...) {
  this$.refTuple;
}, protected=TRUE)





setMethodS3("isPaired", "CopyNumberChromosomalModel", function(this, ...) {
  as.logical(this$.paired);
})




setMethodS3("getReferenceSetTuple", "CopyNumberChromosomalModel", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  cesTuple <- getSetTuple(this);
  cesList <- getListOfSets(cesTuple);

  refTuple <- getRefSetTuple(this);
  if (!force && inherits(refTuple, class(cesTuple)[1])) {
    return(refTuple);
  }

  verbose && enter(verbose, "Building tuple of reference sets");

  # Build a reference set tuple, if missing
  refList <- vector("list", length(cesList));
  for (kk in seq(along=cesList)) {
    ces <- cesList[[kk]];
    refSet <- refList[[kk]];

    if (force || !inherits(refSet, "CopyNumberDataSet")) {
      if (force) {
        verbose && cat(verbose, "Forced recalculation requested.");
      } else {
        verbose && cat(verbose, "No reference available.");
      }
      verbose && enter(verbose, "Calculating average copy-number signals");
      refFile <- getAverageFile(ces, force=force, verbose=less(verbose));
      refSet <- clone(ces);
      clearCache(refSet);
      refFiles <- rep(list(refFile), nbrOfArrays(cesTuple));
      refSet$files <- refFiles;
      rm(refFiles);
      verbose && exit(verbose);
      refList[[kk]] <- refSet;
    }
  } # for (kk ...)

  # Coerce to CopyNumberDataSetTuple, if needed
  refTuple <- as.CopyNumberDataSetTuple(refList);

  this$.referenceTuple <- refTuple;

  verbose && exit(verbose);

  refTuple;
})



setMethodS3("getDataFileMatrix", "CopyNumberChromosomalModel", function(this, array, ..., verbose=FALSE) {
  cesTuple <- getSetTuple(this);
  refTuple <- getReferenceSetTuple(this);

  ceList <- getArrayTuple(cesTuple, array=array, ..., verbose=less(verbose,1));
  rfList <- getArrayTuple(refTuple, array=array, ..., verbose=less(verbose,1));

  # Sanity check
  if (!identical(names(ceList), names(rfList)))
    throw("Internal error. Reference files of non-matching chip types.");

  files <- c(ceList, rfList);
  dim(files) <- c(length(ceList), 2);
  dimnames(files) <- list(names(ceList), c("test", "reference"));
  files;  
}, protected=TRUE);





setMethodS3("getRawCnData", "CopyNumberChromosomalModel", function(this, ceList, refList, chromosome, units=NULL, reorder=TRUE, ..., maxNAFraction=1/5, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ceList':
  if (!inherits(ceList, "list")) {
    throw("Argument 'ceList' is not a list: ", class(ceList)[1]);
  } else {
    for (kk in seq(along=ceList)) {
      ce <- ceList[[kk]];
      className <- "CopyNumberDataFile";
      if (!is.null(ce) && !inherits(ce, className)) {
        throw("Argument 'ceList' contains a non-", className, ": ", 
                                                               class(ce)[1]);
      }
    }
  }

  # Argument 'refList':
  if (!inherits(refList, "list")) {
    throw("Argument 'refList' is not a list: ", class(refList)[1]);
  } else {
    if (length(refList) != length(ceList)) {
      throw("Argument 'refList' is of a different length than 'cesList': ", length(refList), " != ", length(ceList));
    }
    for (kk in seq(along=refList)) {
      ref <- refList[[kk]];
      className <- "CopyNumberDataFile";
      if (!is.null(ref) && !inherits(ref, className)) {
        throw("Argument 'refList' contains a non-", className, ": ", 
                                                              class(ref)[1]);
      }
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving raw CN data");

  # Data set attributes
  chipTypes <- getChipTypes(this);
  arrayNames <- getArrays(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (x, M, stddev, chiptype, unit) from all chip types
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving relative copy-number estimates");
  # Get the chip types as a factor
  chipTypes <- as.factor(chipTypes);
  df <- NULL;
  for (kk in seq(along=chipTypes)) {
    chipType <- chipTypes[kk];
    verbose && enter(verbose, "Chip type: ", chipType);
    ce <- ceList[[kk]];
    if (!is.null(ce)) {
      ref <- refList[[kk]];

      # AD HOC. /HB 2007-09-29
      # Hmmm... what is this? /HB 2009-11-18
      # At least, replaced AffymetrixCelFile 
      # with CopyNumberDataFile. /HB 2009-11-18
      if (!inherits(ref, "CopyNumberDataFile")) {
        ref <- NULL;
      }

      # AD HOC. getXAM() is currently only implemented in the 
      # AffymetrixCelFile class and not defined in any Interface class.
      # It should be implemented by others too, alternatively be 
      # replaced by a "better" method, e.g. extractRawCopyNumbers().
      # /HB 2009-11-18.
      df0 <- getXAM(ce, other=ref, chromosome=chromosome, units=units, 
                                                verbose=less(verbose));
      df0 <- df0[,c("x", "M"), drop=FALSE];
      # Argument 'units' may be NULL
      units0 <- as.integer(rownames(df0));
      verbose && cat(verbose, "Number of units: ", length(units0));

      # Get (theta, sigma) of theta (estimated across all arrays).
      theta <- extractMatrix(ref, units=units0, field="theta", 
                                             drop=TRUE, verbose=verbose);
  
      sdTheta <- extractMatrix(ref, units=units0, field="sdTheta", 
                                             drop=TRUE, verbose=verbose);

      # Estimate the std dev of the raw log2(CN). 
      # [only if ref is averaged across arrays]
      if (isAverageFile(ref)) {
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        # BEGIN: NEED SPECIAL ATTENTION IF ALLELE-SPECIFIC ESTIMATES
        # (whose are not supported yet /HB 2009-11-18)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        # Number of arrays used when averaging (per unit)
        ns <- getNumberOfFilesAveraged(ref, units=units0, verbose=verbose);

        # Sanity check
        stopifnot(length(ns) == length(theta));

        # Use Gauss' approximation (since mu and sigma are on the 
        # intensity scale)
        sdM <- log2(exp(1)) * sqrt(1+1/ns) * sdTheta / theta;

        rm(ns);
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        # END: NEED SPECIAL ATTENTION IF ALLELE-SPECIFIC ESTIMATES
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      } else {
        naValue <- as.double(NA);
        sdM <- rep(naValue, length(theta));
      } # if (isAverageFile(ref))

      verbose && enter(verbose, "Scanning for non-finite values");
      n <- sum(!is.finite(df0[,"M"]));
      fraction <- n / nrow(df0);
      verbose && printf(verbose, "Number of non-finite values: %d (%.1f%%)\n", 
                                                             n, 100*fraction);

      # Sanity check
      if (fraction > maxNAFraction) {
        throw(sprintf("Something is wrong with the data. Too many non-finite values: %d (%.1f%% > %.1f%%)", as.integer(n), 100*fraction, 100*maxNAFraction));
      }
      verbose && exit(verbose);
  
      # Append SD, chip type, and CDF units.
      df0 <- cbind(df0, sdTheta=sdTheta, sdM=sdM, 
                   chipType=rep(chipType, length=length(units0)), unit=units0);
      rm(theta, sdTheta);
  
      df <- rbind(df, df0);
      colnames(df) <- colnames(df0);
      rm(df0, units0);
    } else {
      verbose && cat(verbose, "No copy-number estimates available: ", arrayNames[kk]);
    }

    # Garbage collect
    gc <- gc();
    
    verbose && exit(verbose);
  } # for (kk in ...)
  verbose && exit(verbose);
  

  if (reorder) {
    verbose && enter(verbose, "Re-order by physical position");
    df <- df[order(df[,"x"]),, drop=FALSE];
    rownames(df) <- NULL;
    verbose && exit(verbose);
  }

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  nbrOfUnits <- nrow(df);
  verbose && cat(verbose, sprintf("Extracted data for %d SNPs", nbrOfUnits));
  verbose && exit(verbose);

  df;
}, protected=TRUE)



setMethodS3("calculateChromosomeStatistics", "CopyNumberChromosomalModel", function(this, arrays=NULL, chromosomes=getChromosomes(this), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  allChromosomes <- getChromosomes(this);

  # Argument 'arrays':
  arrays <- indexOfArrays(this, arrays=arrays);

  # Argument 'chromosomes':
  if (is.null(chromosomes)) {
    chromosomes <- allChromosomes;
  } else {
    unknown <- chromosomes[!(chromosomes %in% allChromosomes)];
    if (length(unknown) > 0) {
      throw("Argument 'chromosomes' contains unknown values: ", paste(unknown, collapse=", "));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimating number of copies per chromosome");
  nbrOfChromosomes <- length(chromosomes);
  verbose && printf(verbose, "Chromosomes (%d): %s", nbrOfChromosomes, paste(chromosomes, collapse=", "));

  res <- list();
  arrayNames <- getNames(this)[arrays];
  nbrOfArrays <- length(arrayNames);

  for (aa in seq(length=nbrOfArrays)) {
    array <- arrays[aa];
    arrayName <- arrayNames[aa];

    files <- getDataFileMatrix(this, array=array, verbose=less(verbose,5));
    ceList <- files[,"test"];
    rfList <- files[,"reference"];

    res[[arrayName]] <- list();
    for (chromosome in chromosomes) {
      verbose && enter(verbose, 
                         sprintf("Array #%d ('%s') of %d on chromosome %s", 
                                  aa, arrayName, nbrOfArrays, chromosome));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Get (x, M, stddev, chiptype, unit) data from all chip types
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      data <- getRawCnData(this, ceList=ceList, refList=rfList, 
                        chromosome=chromosome, ..., verbose=less(verbose));
      M <- data[,"M"];
      rm(data);

      fit <- list(
        mean = mean(M, na.rm=TRUE),
        sd = sd(M, na.rm=TRUE),
        median = median(M, na.rm=TRUE),
        mad = mad(M, na.rm=TRUE),
        quantiles = quantile(M, probs=seq(0,1,by=0.01), na.rm=TRUE),
        sdDiff = sd(diff(M), na.rm=TRUE)/sqrt(2),
        madDiff = mad(diff(M), na.rm=TRUE)/sqrt(2),
        nbrOfLoci = length(M),
        nbrOfNAs = sum(is.na(M))
      );

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Estimate
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Garbage collection
      gc <- gc();
      verbose && print(verbose, gc);

      res[[arrayName]][[chromosome]] <- fit;
      verbose && exit(verbose);
    } # for (chromosome in ...)
  } # for (aa in ...)
  verbose && exit(verbose);

  res;
}, protected=TRUE)  # calculateChromosomeStatistics()




###########################################################################/**
# @RdocMethod fit
#
# @title "Fits the model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{See subclasses.}
# }
#
# \value{
#  See subclasses.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fit", "CopyNumberChromosomalModel", abstract=TRUE);



setMethodS3("extractRawCopyNumbers", "CopyNumberChromosomalModel", function(this, array, chromosome, ..., cache=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'array':
  array <- Arguments$getIndex(array, range=c(1,nbrOfArrays(this)));

  # Argument 'chromosome':
  allChromosomes <- getChromosomes(this);
  chromosome <- Arguments$getIndex(chromosome, range=range(allChromosomes));
  if (!chromosome %in% allChromosomes)
    throw("Argument 'chromosome' has an unknown value: ", chromosome);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  key <- list(method="extractRawCopyNumbers", class=class(this), array=array, chromosome=chromosome, ...);
  id <- digest(key);
  cacheList <- this$.extractRawCopyNumbersCache;
  if (!is.list(cacheList))
    cacheList <- list();
  rawCNs <- cacheList[[id]];
  if (!force && !is.null(rawCNs)) {
    return(rawCNs);
  }

  # Extract the test and reference arrays
  files <- getDataFileMatrix(this, array=array, verbose=less(verbose,5));
  ceList <- files[,"test"];
  rfList <- files[,"reference"];

  data <- getRawCnData(this, ceList=ceList, refList=rfList, 
                          chromosome=chromosome, ..., verbose=less(verbose));

  rawCNs <- RawCopyNumbers(cn=data[,"M"], x=data[,"x"], 
                           chromosome=chromosome); 

  # Save to cache?
  if (cache) {
    cacheList[[id]] <- rawCNs;
    this$.extractRawCopyNumbersCache <- cacheList;
    rm(cacheList);
  }

  rawCNs;
})




###########################################################################/**
# @RdocMethod estimateSds
#
# @title "Estimates the standard deviation of the raw copy numbers (log2-ratios) robustly"
#
# \description{
#  @get "title" using a first-order difference variance estimator, which is
#  an estimator that is fairly robust for change points.
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{The arrays to be queried.}
#   \item{chromosomes}{The chromosomes to be queried.}
#   \item{...}{Additional arguments passed to 
#      @seemethod "extractRawCopyNumbers".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a CxK @double @matrix where C is the number of chromosomes, 
#  and K is the number of arrays (arrays are always the last dimension).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/########################################################################### 
setMethodS3("estimateSds", "CopyNumberChromosomalModel", function(this, arrays=seq(length=nbrOfArrays(this)), chromosomes=getChromosomes(this), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  arrays <- indexOfArrays(this, arrays=arrays);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  nbrOfChromosomes <- length(chromosomes);

  naValue <- as.double(NA);
  res <- matrix(naValue, nrow=nbrOfChromosomes, ncol=length(arrays));
  colnames(res) <- getArrays(this)[arrays];
  rownames(res) <- chromosomes;


  for (rr in seq(length=nbrOfChromosomes)) {
    chromosome <- chromosomes[rr];
    verbose && enter(verbose, "Chromosome #", rr, " of ", nbrOfChromosomes);

    for (cc in seq(along=arrays)) {
      array <- arrays[cc];
      verbose && enter(verbose, "Array #", cc, " of ", length(arrays));

      rawCns <- extractRawCopyNumbers(this, array=array, chromosome=chromosome, ..., verbose=less(verbose,5));

#      verbose && enter(verbose, "First-order robust variance estimator");
      res[rr,cc] <- estimateStandardDeviation(rawCns);
      rm(rawCns);
#      verbose && exit(verbose);

      verbose && exit(verbose);
    }

    verbose && exit(verbose);
  }

  res;
}, protected=TRUE)




##############################################################################
# HISTORY:
# 2009-11-18
# o CLEAN UP: Removed all Affymetrix specific classes/methods; using Interface 
#   classes almost everywhere.  It is a bit ad hoc design, but we will 
#   worry about that later; let's get something working first.
# o BUG FIX: getRawCnData(..., reorder=FALSE, verbose=TRUE) would give an
#   error because 'nbrOfUnits' would not have  been defined.
# 2009-11-16
# o CLEAN UP: Using getDataFileMatrix() instead of the old name
#   getMatrixChipEffectFiles().
# 2009-11-13
# o ROBUSTNESS: Now arguments 'ces' and 'ref' and CopyNumberChromosomalModel
#   have to be CnChipEffectFile|Set, otherwise an exception is thrown.
#   Before it was possible to pass a SnpChipEffectSet unnoticed, although
#   only total CNs are modelled.  Thanks Pierre Neuvial for this report.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: estimateSds().
# 2008-07-16
# o Added support for specifying the reference by a single file (or a list of
#   files if more than one set is used).
# 2008-07-01
# o MEMORY OPTIMIZATION: When calling extractRawCopyNumbers(obj) on an
#   CopyNumberChromosomalModel object, the result would be cached in memory
#   (in the object). This would result in an increasing memory usage when
#   data was extracted from more and more arrays. The cache could be cleared
#   by calling gc(obj), but avoid this problem by default, the method does
#   no longer cache results.  To cache, the method has to be called with
#   argument 'cache=TRUE'.  Thanks Jason Li for reporting this.
# 2008-06-07
# o BUG FIX: getReferenceSetTuple() of CopyNumberChromosomalModel would
#   generate nbrOfArrays(ces) instead of nbrOfArrays(cesTuple) reference
#   files if not reference set was specified.
# 2008-05-31
# o BUG FIX: getRawCnData() of CopyNumberChromosomalModel threw "Exception: 
#   Argument 'ceList' contains a non-ChipEffectFile: NULL" if multiple
#   ChipEffectSet:s is modelled and one of them don't have all arrays.
#   Thanks Lavinia Gordon for spotting this.
# 2008-05-08
# o BUG FIX: getRawCnData() of CopyNumberChromosomalModel gave an error if
# 2008-03-10
# o Added estimateSds() with Rdoc comments.
# 2007-11-27
# o Changed default 'maxNAFraction' to 1/5 (from 1/8) in getRawCnData().
# o BUG FIX: Two different clearCache() was defined.
# 2007-10-17
# o Added extractRawCopyNumbers().
# o Created.
##############################################################################
