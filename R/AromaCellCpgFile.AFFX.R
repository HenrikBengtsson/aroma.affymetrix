
setMethodS3("allocateFromCdf", "AromaCellCpgFile", function(static, cdf, path=getPath(cdf), tags="*", ...) {
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
})



############################################################################
# HISTORY:
# 2010-02-19 [MR]
# o Created from AromaCellPositionFile.AFFX.R.
############################################################################
