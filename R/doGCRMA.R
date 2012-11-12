#  doGCRMA() runs in bounded memory and replicates the results of
#  @see "gcrma::gcrma" in the \pkg{gcrma} package with great precision.
# \references{
#  [1] Z. Wu, R. Irizarry, R. Gentleman, F.M. Murillo & F. Spencer, A Model Based Background Adjustment for Oligonucleotide Expression Arrays, JASA, 2004.
# }

setMethodS3("doGCRMA", "AffymetrixCelSet", function(csR, arrays=NULL, type=c("fullmodel", "affinities"), ..., uniquePlm=FALSE, drop=TRUE, ram=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'csR':
  className <- "AffymetrixCelSet";
  if (!inherits(csR, className)) {
    throw(sprintf("Argument 'csR' is not a %s: %s", className, class(csR)[1]));
  }

  # Argument 'arrays':
  if (!is.null(arrays)) {
    throw("Not supported. Argument 'arrays' should be NULL.");
    arrays <- Arguments$getIndices(arrays, max=length(csR));
  }

  # Argument 'uniquePlm':
  uniquePlm <- Arguments$getLogical(uniquePlm);

  # Argument 'drop':
  drop <- Arguments$getLogical(drop);

  # Argument 'type':
  type <- match.arg(type);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);




  verbose && enter(verbose, "GCRMA");
  verbose && cat(verbose, "Arguments:");
  arraysTag <- seqToHumanReadable(arrays);
  verbose && cat(verbose, "arrays:");
  verbose && str(verbose, arraysTag);
  verbose && cat(verbose, "Fit PLM on unique probe sets: ", uniquePlm);
  verbose && cat(verbose, "ram: ", ram);


  # List of objects to be returned
  res <- list();
  if (!drop) {
    res <- c(res, list(csR=csR));
  }

  verbose && cat(verbose, "Data set");
  verbose && print(verbose, csR);

  if (!is.null(arrays)) {
    verbose && enter(verbose, "GCRMA/Extracting subset of arrays");
    csR <- extract(csR, arrays);
    verbose && cat(verbose, "Data subset");
    verbose && print(verbose, csR);
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "GCRMA/Background correction");

  # Currently, you must use the standard CDF file. 
  cdf <- getCdf(csR);
  chipTypeS <- getChipType(cdf, fullname=FALSE);
  cdfS <- AffymetrixCdfFile$byChipType(chipTypeS);
  if (!equals(cdfS, cdf)) {
    setCdf(csR, cdfS);
    on.exit({
      # Make sure to undo, e.g. if interrupted.
      setCdf(csR, cdf);
    });
  }

  bc <- GcRmaBackgroundCorrection(csR, type=type);
  verbose && print(verbose, bc);
  csB <- process(bc, verbose=verbose);
  verbose && print(verbose, csB);

  if (!equals(cdfS, cdf)) {
    # Now, use the custom CDF in what follows
    setCdf(csB, cdf);
  }
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(bc=bc, csB=csB));
  }

  # Clean up
  rm(csR, bc);
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "GCRMA/Rank-based quantile normalization (PM-only)");
  qn <- QuantileNormalization(csB, typesToUpdate="pm");
  verbose && print(verbose, qn);
  csN <- process(qn, verbose=verbose);
  verbose && print(verbose, csN);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(qn=qn, csN=csN));
  }

  # Clean up
  rm(csB, qn);
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "GCRMA/Probe summarization (log-additive model fitted using median polish)");
  verbose && cat(verbose, "Fit PLM on unique probe sets: ", uniquePlm);

  if (uniquePlm) {
    verbose && enter(verbose, "Probe-summarization using a \"unique\" CDF requested");

    verbose && enter(verbose, "Getting \"unique\" CDF (with non-unique probes dropped)")
    cdf <- getCdf(csN);
    verbose && cat(verbose, "CDF:");
    verbose && print(verbose, cdf);
    cdfU <- getUniqueCdf(cdf, verbose=less(verbose, 5));
    verbose && cat(verbose, "CDF with non-unique probes dropped:");
    verbose && print(verbose, cdfU);
    verbose && exit(verbose)

    if (equals(cdfU, cdf)) {
      verbose && cat(verbose, "The \"unique\" CDF equals the original CDF: Skipping.");
    } else {
      csNU <- convertToUnique(csN, verbose=verbose);
      verbose && print(verbose, csNU);
      csN <- csNU;
    }
    verbose && exit(verbose);
  }

  plm <- RmaPlm(csN, flavor="oligo");
  verbose && print(verbose, plm);
  if (length(findUnitsTodo(plm)) > 0) {
    units <- fit(plm, ram=ram, verbose=verbose);
    verbose && str(verbose, units);
    rm(units);
  }
  verbose && print(verbose, gc);
  ces <- getChipEffectSet(plm);
  verbose && print(verbose, ces);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(ces=ces, plm=plm));
  }

  # Clean up
  rm(plm, csN);
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  # Return only the final output data set?
  if (drop) {
    res <- ces;
  }

  res;
}) # doGCRMA()


setMethodS3("doGCRMA", "character", function(dataSet, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "GCRMA");

  verbose && enter(verbose, "GCRMA/Setting up CEL set");
  csR <- AffymetrixCelSet$byName(dataSet, ..., verbose=less(verbose, 50),
                                                  .onUnknownArgs="ignore");
  verbose && print(verbose, csR);
  verbose && exit(verbose);

  res <- doGCRMA(csR, ..., verbose=verbose);

  # Clean up
  rm(csR);
  gc <- gc();

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2011-11-10
# o ROBUSTNESS: doGCRMA() is now guaranteed to undo any changes of 
#   the CDF of the data set, e.g. if there is a user interrupt.
# 2011-04-07
# o Added argument 'drop'.
# o Added argument 'uniquePlm'.
# 2010-09-26
# o Now doGCRMA() automagically makes sure that the default CDF is used
#   in the GcRmaBackgroundCorrection step, while use a custom CDF
#   everywhere else if set.
# o Added argument 'type' to doGCRMA() which is passed to 
#   QuantileNormalization().
# 2010-08-14
# o Created from doRMA.R.
############################################################################
