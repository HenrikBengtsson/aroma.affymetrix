# \section{Reproducibility of affy}{
#   This implementation of the RMA method reproduces @see "affy::justRMA"
#   in \pkg{affy} package quite well.  It does so by still using a
#   constant memory profile, i.e. it is possible to use this implementation
#   to run RMA on a much large data set than what is possible with
#   \pkg{affy}.  At least 20-50 \emph{times} more samples should be
#   doable, if not more.
# }
setMethodS3("justRMA", "AffymetrixCelSet", function(csR, flavor=c("oligo", "affyPLM"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csR':
  className <- "AffymetrixCelSet";
  if (!inherits(csR, className)) {
    throw(sprintf("Argument 'csR' is not a %s: %s", className, class(csR)[1]));
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "justRMA() via aroma.affymetrix");

  # Run RMA
  ces <- doRMA(csR, drop=TRUE, flavor=flavor, ..., verbose=verbose);

  # Extract estimates
  eset <- extractExpressionSet(ces, orderUnitsBy="lexicographic",
                               annotationPkg="ChipDb", verbose=verbose);

  # BACKWARD COMPATIBIILTY/AD HOC:
  # Coerce std errors to missing values in case they are all 1.
  se.exprs <- assayData(eset)$se.exprs
  if (all(se.exprs == 1)) {
    is.na(se.exprs) <- TRUE;
    assayData(eset)$se.exprs <- se.exprs;
  }
  se.exprs <- NULL;

  # Set the protocol data; scan dates
  data <- data.frame(ScanDate=getTimestamps(csR));
  rownames(data) <- gsub(",chipEffects", "", getFullNames(ces), fixed=TRUE);
  data <- as(data, "AnnotatedDataFrame");
  protocolData(eset) <- data;
  data <- NULL;

  verbose && print(verbose, eset);
  verbose && exit(verbose);

  eset;
}, protected=TRUE) # justRMA()

# HISTORY:
# 2014-04-27
# o Added  ...finally!
