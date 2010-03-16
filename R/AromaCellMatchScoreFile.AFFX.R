setMethodS3("allocateFromCdf", "AromaCellMatchScoreFile", function(static, cdf, path=getPath(cdf), tags="*", ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=NULL);

  chipType <- getChipType(cdf);
  parts <- strsplit(chipType, split=",", fixed=TRUE);
  parts <- unlist(parts, use.names=FALSE);

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
# 2010-03-14 [HB]
# o BUG FIX: allocateFromCdf() of AromaCellMatchScoreFile would drop all
#   but the first tag.
# 2009-02-22 [HB]
# o Forgot to make allocateFromCdf() of AromaCellMatchScoreFile static.
# 2009-02-10 [HB]
# o Added optional validation of number of cells to byChipType().
# o Static method byChipType() was not declared static.
# 2008-10-28 [MR]
# o Created from AromaCellSequenceFile.R.
############################################################################
