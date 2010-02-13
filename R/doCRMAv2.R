setMethodS3("doCRMAv2", "AffymetrixCelSet", function(csR, combineAlleles=TRUE, arrays=NULL, ..., verbose=FALSE) {
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

  # Argument 'arrays':
  if (!is.null(arrays)) {
    arrays <- Arguments$getIndices(arrays, max=nbrOfArrays(csR));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "CRMAv2");
  verbose && cat(verbose, "Arguments:");
  verbose && cat(verbose, "combineAlleles: ", combineAlleles);
  verbose && cat(verbose, "arrays: ", arrays);

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
  
  verbose && enter(verbose, "CRMAv2/Base position normalization");
  bpn <- BasePositionNormalization(csC, target="zero");
  verbose && print(verbose, bpn);
  csN <- process(bpn, verbose=verbose);
  verbose && print(verbose, csN);
  verbose && exit(verbose);
  
  verbose && enter(verbose, "CRMAv2/Probe summarization");
  plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=combineAlleles);
  verbose && print(verbose, plm);
  if (length(findUnitsTodo(plm)) > 0) {
    # Fit CN probes quickly (~5-10s/array + some overhead)
    units <- fitCnProbes(plm, verbose=verbose);
    verbose && str(verbose, units);
    # Fit remaining units, i.e. SNPs (~5-10min/array)
    units <- fit(plm, verbose=verbose);
    verbose && str(verbose, units);
  }  
  ces <- getChipEffectSet(plm);
  verbose && print(verbose, ces);
  verbose && exit(verbose);
  
  verbose && enter(verbose, "CRMAv2/PCR fragment-length normalization");
  fln <- FragmentLengthNormalization(ces, target="zero");
  verbose && print(verbose, fln);
  cesN <- process(fln, verbose=verbose);
  verbose && print(verbose, cesN);
  verbose && exit(verbose);
  
  verbose && enter(verbose, "CRMAv2/Export to technology-independent data files");
  dsN <- exportTotalAndFracB(cesN, verbose=verbose);
  verbose && print(verbose, dsN);
  verbose && exit(verbose);

  verbose && exit(verbose);

  dsN;
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

  dsN <- doCRMAv2(csR, ..., verbose=verbose);

  verbose && exit(verbose);

  dsN;
})


############################################################################
# HISTORY:
# 2010-02-13
# o Restructured.  Now there is a doCRMAv2() for AffymetrixCelSet:s and
#   one for character strings.
# 2009-09-13
# o (Re)created.
############################################################################