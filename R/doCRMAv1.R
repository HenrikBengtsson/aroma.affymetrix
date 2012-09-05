setMethodS3("doCRMAv1", "AffymetrixCelSet", function(csR, shift=+300, combineAlleles=TRUE, lengthRange=NULL, arrays=NULL, ..., drop=TRUE, ram=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'csR':
  className <- "AffymetrixCelSet";
  if (!inherits(csR, className)) {
    throw(sprintf("Argument 'csR' is not a %s: %s", className, class(csR)[1]));
  }

  # Argument 'shift':
  shift <- Arguments$getNumeric(shift);

  # Argument 'combineAlleles':
  combineAlleles <- Arguments$getLogical(combineAlleles);

  # Argument 'arrays':
  if (!is.null(arrays)) {
    throw("Not supported. Argument 'arrays' should be NULL.");
    arrays <- Arguments$getIndices(arrays, max=nbrOfArrays(csR));
  }

  # Argument 'drop':
  drop <- Arguments$getLogical(drop);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "CRMAv1");
  verbose && cat(verbose, "Arguments:");
  verbose && cat(verbose, "combineAlleles: ", combineAlleles);
  arraysTag <- seqToHumanReadable(arrays);
  verbose && cat(verbose, "arrays:");
  verbose && str(verbose, arraysTag);
  verbose && cat(verbose, "ram: ", ram);


  # List of objects to be returned
  res <- list();
  if (!drop) {
    res <- c(res, list(csR=csR));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Setup data set to be processed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && cat(verbose, "Data set");
  verbose && print(verbose, csR);

  if (!is.null(arrays)) {
    verbose && enter(verbose, "CRMAv1/Extracting subset of arrays");
    csR <- extract(csR, arrays);
    verbose && cat(verbose, "Data subset");
    verbose && print(verbose, csR);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # CRMAv1
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "CRMAv1/Allelic crosstalk calibration");
  acc <- AllelicCrosstalkCalibration(csR, model="CRMA", tags="*,v1");
  verbose && print(verbose, acc);
  csC <- process(acc, verbose=verbose);
  verbose && print(verbose, csC);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(acc=acc, csC=csC));
  }

  # Clean up
  rm(csR, acc);
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "CRMAv1/Probe summarization");
  plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=combineAlleles, 
                                                            shift=shift);
  verbose && print(verbose, plm);
  if (length(findUnitsTodo(plm)) > 0) {
    # Fit CN probes quickly (~5-10s/array + some overhead)
    units <- fitCnProbes(plm, verbose=verbose);
    verbose && str(verbose, units);
    # Fit remaining units, i.e. SNPs (~5-10min/array)
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
  rm(plm, csC);
  gc <- gc();
  
  verbose && enter(verbose, "CRMAv1/PCR fragment-length normalization");
  fln <- FragmentLengthNormalization(ces, target="zero", lengthRange=lengthRange);
  verbose && print(verbose, fln);
  cesN <- process(fln, verbose=verbose);
  verbose && print(verbose, cesN);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(fln=fln, cesN=cesN));
  }

  # Clean up
  rm(fln, ces);
  gc <- gc();
  
  verbose && enter(verbose, "CRMAv1/Export to technology-independent data files");
  dsNList <- exportTotalAndFracB(cesN, verbose=verbose);
  verbose && print(verbose, dsNList);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(dsNList=dsNList));
  }

  # Clean up
  rm(cesN);
  gc <- gc();

  verbose && exit(verbose);

  # Return only the final results?
  if (drop) {
    res <- dsNList;
  }

  res;
}) # doCRMAv1()


setMethodS3("doCRMAv1", "character", function(dataSet, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "CRMAv1");

  verbose && enter(verbose, "CRMAv1/Setting up CEL set");
  csR <- AffymetrixCelSet$byName(dataSet, ..., verbose=less(verbose, 50),
                                                  .onUnknownArgs="ignore");
  verbose && print(verbose, csR);
  verbose && exit(verbose);

  dsNList <- doCRMAv1(csR, ..., verbose=verbose);

  # Clean up
  rm(csR);
  gc <- gc();

  verbose && exit(verbose);

  dsNList;
})


setMethodS3("doASCRMAv1", "default", function(...) {
  doCRMAv1(..., combineAlleles=FALSE); 
})


############################################################################
# HISTORY:
# 2012-09-05
# o ROBUSTNESS: Now doCRMAv1() adds also tag "v1" to the allele-specific
#   calibration step.  The reason for this is to differentiate it from 
#   the output of doCRMAv2().  NOTE: This update means that any old CRMAv1
#   analyses will not be detected by doCRMAv1(); to have doCRMAv1() detect 
#   those add tag "v1" in that calibration step, e.g. "ACC,-XY,v1".
# 2011-04-07
# o Added argument 'drop'.
# 2011-03-14
# o doCRMAv1() gained argument 'lengthRange', which is passed to 
#   the constructor of FragmentLengthNormalization.
# 2010-06-21
# o Added doASCRMAv1() for a convenient allele-specific CRMAv1 wrapper.
# 2010-06-07
# o Added argument shift=+300 to doCRMAv1().
# 2010-05-17
# o CORRECTION: doCRMAv1() forgot to shift +300 the signals before 
#   doing the probe-level summarization.
# 2010-04-21
# o BUG FIX: doCRMAv1() for AffymetrixCelSet used undefined 'csN' internally
#   instead of 'csC'.
# 2010-04-04
# o Created from doCRMAv2.R.
# o (Re)created.
############################################################################
