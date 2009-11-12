setMethodS3("extractSnpFeatureSet", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  pathnames <- getPathnames(this);

  verbose && enter(verbose, "Reading complete SnpFeatureSet");
  verbose && cat(verbose, "Number of data files: ", length(this));
  verbose && cat(verbose, "Chip type: ", getChipType(this));
  verbose2 <- as.logical(verbose);
  res <- read.celfiles(pathnames, ..., verbose=verbose2);
  verbose && exit(verbose);

  res;
}) # extractSnpFeatureSet()


############################################################################
# HISTORY:
# 2009-10-16
# o Created.
############################################################################ 
