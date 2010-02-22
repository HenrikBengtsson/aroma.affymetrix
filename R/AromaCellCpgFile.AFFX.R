setMethodS3("allocateFromCdf", "AromaCellCpgFile", function(static, cdf, path=getPath(cdf), tags="*", ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  tags <- strsplit(tags, split=",", fixed=TRUE)[[1]];
  chipType <- getChipType(cdf);
  parts <- strsplit(chipType, split=",", fixed=TRUE)[[1]];
  tags[tags == "*"] <- paste(parts[-1], collapse=",");
  chipType <- parts[1];

  # Drop unwanted tags. AD HOC. /HB 2010-02-22
  chipType <- gsub(",monocell", "", chipType, fixed=TRUE);

  platform <- getPlatform(cdf);
  nbrOfCells <- nbrOfCells(cdf);
  fullname <- paste(c(chipType, tags), collapse=",");
  ext <- getFilenameExtension(static);
  filename <- sprintf("%s.%s", fullname, ext);
  allocate(static, filename=filename, path=path, nbrOfCells=nbrOfCells, 
                                 platform=platform, chipType=chipType, ...);
}, static=TRUE)



############################################################################
# HISTORY:
# 2009-02-22 [HB]
# o Forgot to make allocateFromCdf() of AromaCellCpgFile static.
# 2010-02-19 [MR]
# o Created from AromaCellPositionFile.AFFX.R.
############################################################################
