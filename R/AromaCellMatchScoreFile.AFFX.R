# @author "MR, HB"
setMethodS3("allocateFromCdf", "AromaCellMatchScoreFile", function(static, cdf, path=NULL, tags="*", ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile")

  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=NULL)

  chipType <- getChipType(cdf)
  parts <- strsplit(chipType, split=",", fixed=TRUE)
  parts <- unlist(parts, use.names=FALSE)

  tags[tags == "*"] <- paste(parts[-1], collapse=",")
  chipType <- parts[1]
  chipType <- gsub(",monocell", "", chipType)

  # Output path
  if (is.null(path)) {
    path <- file.path("annotationData", "chipTypes", chipType)
  }
  path <- Arguments$getWritablePath(path)

  platform <- getPlatform(cdf)
  nbrOfCells <- nbrOfCells(cdf)
  fullname <- paste(c(chipType, tags), collapse=",")
  ext <- getFilenameExtension(static)
  filename <- sprintf("%s.%s", fullname, ext)
  allocate(static, filename=filename, path=path, nbrOfCells=nbrOfCells,
                   platform=platform, chipType=chipType, ...)
}, static=TRUE)
