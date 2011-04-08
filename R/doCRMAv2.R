setMethodS3("doCRMAv2", "AffymetrixCelSet", function(csR, combineAlleles=TRUE, lengthRange=NULL, arrays=NULL, ..., plm=c("AvgCnPlm", "RmaCnPlm"), drop=TRUE, ram=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'csR':
  className <- "AffymetrixCelSet";
  if (!inherits(csR, className)) {
    throw(sprintf("Argument 'csR' is not a %s: %s", className, class(csR)[1]));
  }

  # Argument 'combineAlleles':
  combineAlleles <- Arguments$getLogical(combineAlleles);

  # Argument 'plm':
  plm <- match.arg(plm);

  # Argument 'arrays':
  if (!is.null(arrays)) {
    arrays <- Arguments$getIndices(arrays, max=nbrOfArrays(csR));
    if (plm == "RmaCnPlm") {
      throw(sprintf("Argument 'arrays' must not be specified when argument 'plm' is \"%s\".", plm));
    }
  }

  # Argument 'drop':
  drop <- Arguments$getLogical(drop);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "CRMAv2");
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


  verbose && cat(verbose, "Data set");
  verbose && print(verbose, csR);

  if (!is.null(arrays)) {
    verbose && enter(verbose, "CRMAv2/Extracting subset of arrays");
    csR <- extract(csR, arrays);
    verbose && cat(verbose, "Data subset");
    verbose && print(verbose, csR);
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "CRMAv2/Allelic crosstalk calibration");
  acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2");
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

  verbose && enter(verbose, "CRMAv2/Base position normalization");
  bpn <- BasePositionNormalization(csC, target="zero");
  verbose && print(verbose, bpn);
  csN <- process(bpn, verbose=verbose);
  verbose && print(verbose, csN);
  verbose && exit(verbose);

  # Clean up
  rm(csC, bpn);
  gc <- gc();
  verbose && print(verbose, gc);
  
  verbose && enter(verbose, "CRMAv2/Probe summarization");
  cnPlm <- AvgCnPlm;
  if (plm == "RmaCnPlm") {
    cnPlm <- RmaCnPlm;
  }
  plm <- cnPlm(csN, mergeStrands=TRUE, combineAlleles=combineAlleles);
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
  rm(plm, csN);
  gc <- gc();
  
  verbose && enter(verbose, "CRMAv2/PCR fragment-length normalization");
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
  
  verbose && enter(verbose, "CRMAv2/Export to technology-independent data files");
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
}) # doCRMAv2()


setMethodS3("doCRMAv2", "character", function(dataSet, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "CRMAv2");

  verbose && enter(verbose, "CRMAv2/Setting up CEL set");
  csR <- AffymetrixCelSet$byName(dataSet, ..., verbose=less(verbose, 50),
                                                  .onUnknownArgs="ignore");
  verbose && print(verbose, csR);
  verbose && exit(verbose);

  dsNList <- doCRMAv2(csR, ..., verbose=verbose);

  # Clean up
  rm(csR);
  gc <- gc();

  verbose && exit(verbose);

  dsNList;
})


setMethodS3("doASCRMAv2", "default", function(...) {
  doCRMAv2(..., combineAlleles=FALSE); 
})


############################################################################
# HISTORY:
# 2011-04-07
# o Added argument 'drop'.
# 2011-03-14
# o doCRMAv2() gained argument 'lengthRange', which is passed to 
#   the constructor of FragmentLengthNormalization.
# 2010-06-21
# o Added doASCRMAv2() for a convenient allele-specific CRMAv2 wrapper.
# 2010-04-04
# o Added argument 'plm' to doCRMAv2().
# 2010-02-15
# o MEMORY OPTIMIZATION: Now doCRMAv2() removes as much as possible.
# 2010-02-13
# o Restructured.  Now there is a doCRMAv2() for AffymetrixCelSet:s and
#   one for character strings.
# 2009-09-13
# o (Re)created.
############################################################################
