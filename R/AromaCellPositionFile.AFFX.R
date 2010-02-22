setMethodS3("allocateFromCdf", "AromaCellPositionFile", function(static, cdf, path=getPath(cdf), tags="*", ...) {
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  tags <- strsplit(tags, split=",", fixed=TRUE)[[1]];
  chipType <- getChipType(cdf);
  parts <- strsplit(chipType, split=",", fixed=TRUE)[[1]];
  tags[tags == "*"] <- paste(parts[-1], collapse=",");
  chipType <- parts[1];
  chipType <- gsub(",monocell", "", chipType);
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
# o Forgot to make allocateFromCdf() of AromaCellPositionFile static.
# 2009-02-16 [HB]
# Removed argument 'validate' from byChipType() of AromaCellPositionFile.
# 2009-02-10 [HB]
# o Added optional validation of number of cells to byChipType().
# o Static method byChipType() was not declared static.
# 2008-12-09 [MR]
# o Created from AromaCellMatchScoresFile.R.
############################################################################
