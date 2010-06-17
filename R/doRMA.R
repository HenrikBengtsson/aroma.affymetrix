#  doRMA() runs in bounded memory and replicates the results of
#  fitPLM() in the affyPLM package with great precision.

setMethodS3("doRMA", "AffymetrixCelSet", function(csR, arrays=NULL, ..., ram=NULL, verbose=FALSE) {
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
    arrays <- Arguments$getIndices(arrays, max=nbrOfArrays(csR));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "RMA");
  verbose && cat(verbose, "Arguments:");
  arraysTag <- seqToHumanReadable(arrays);
  verbose && cat(verbose, "arrays:");
  verbose && str(verbose, arraysTag);
  verbose && cat(verbose, "ram: ", ram);


  verbose && cat(verbose, "Data set");
  verbose && print(verbose, csR);

  if (!is.null(arrays)) {
    verbose && enter(verbose, "RMA/Extracting subset of arrays");
    csR <- extract(csR, arrays);
    verbose && cat(verbose, "Data subset");
    verbose && print(verbose, csR);
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "RMA/Background correction (normal & exponential mixture model)");
  bc <- RmaBackgroundCorrection(csR);
  verbose && print(verbose, bc);
  csB <- process(bc, verbose=verbose);
  verbose && print(verbose, csB);
  verbose && exit(verbose);

  # Clean up
  rm(csR, bc);
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "RMA/Rank-based quantile normalization (PM-only)");
  qn <- QuantileNormalization(csB, typesToUpdate="pm");
  verbose && print(verbose, qn);
  csN <- process(qn, verbose=verbose);
  verbose && print(verbose, csN);
  verbose && exit(verbose);

  # Clean up
  rm(csB, qn);
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "RMA/Probe summarization (log-additive model)");
  plm <- RmaPlm(csN);
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

  # Clean up
  rm(plm, csN);
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  ces;
}) # doRMA()


setMethodS3("doRMA", "character", function(dataSet, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "RMA");

  verbose && enter(verbose, "RMA/Setting up CEL set");
  csR <- AffymetrixCelSet$byName(dataSet, ..., verbose=less(verbose, 50),
                                                  .onUnknownArgs="ignore");
  verbose && print(verbose, csR);
  verbose && exit(verbose);

  res <- doRMA(csR, ..., verbose=verbose);

  # Clean up
  rm(csR);
  gc <- gc();

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2010-06-16
# o Created from doCRMAv1.R.
############################################################################
