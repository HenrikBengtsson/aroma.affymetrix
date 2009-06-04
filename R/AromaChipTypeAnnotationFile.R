###########################################################################/**
# @RdocClass AromaChipTypeAnnotationFile
#
# @title "The AromaChipTypeAnnotationFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaChipTypeAnnotationFile object represents an annotation file for a
#  specific chip type.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AromaChipTypeAnnotationFile", function(...) {
  this <- extend(AffymetrixFile(...), "AromaChipTypeAnnotationFile");

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})


setMethodS3("as.character", "AromaChipTypeAnnotationFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Filename: %s", getFilename(this)));
  s <- c(s, sprintf("Filesize: %.2fMB", getFileSize(this)/1024^2));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)



###########################################################################/**
# @RdocMethod fromFile
#
# @title "Sets up an AromaChipTypeAnnotationFile"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filename}{The filename of to the file.}
#  \item{path}{The path to the file.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns an instance of @see "AromaChipTypeAnnotationFile" or its subclasses.
#  If the file is not found or if it is of the wrong file format, an
#  error is thrown.
# }
#
# @author
#
# \seealso{
#   @seemethod "byChipType".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("fromFile", "AromaChipTypeAnnotationFile", function(static, filename, path=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Arguments 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path, 
                                                              mustExist=TRUE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Try to define an instance of a subclass traversing bottom up.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  clazz <- Class$forName(class(static)[1]);
  for (className in rev(getKnownSubclasses(clazz))) {
    clazz <- Class$forName(className);
    tryCatch({
      res <- newInstance(clazz, pathname);
      return(res);
    }, error = function(ex) {})
  }

  newInstance(static, pathname);
}, static=TRUE, protected=TRUE)




###########################################################################/**
# @RdocMethod byChipType
# @aliasmethod byName
#
# @title "Defines an AromaChipTypeAnnotationFile object by chip type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chipType}{A @character string.}
#  \item{tags}{An optional @character @vector of tags.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AromaChipTypeAnnotationFile" object.  
# }
#
# @author
#
# \seealso{
#   @seemethod "fromFile".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("byChipType", "AromaChipTypeAnnotationFile", function(static, chipType, tags=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search for a matching annotation file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- findByChipType(static, chipType, tags=tags);
  if (is.null(pathname)) {
    throw("Could not create ", class(static)[1], " object. No annotation chip type file with that chip type found: ", chipType);
  }

  res <- fromFile(static, filename=pathname, path=NULL, ...);
  verbose && print(verbose, res);

  res;
}, static=TRUE)


setMethodS3("byName", "AromaChipTypeAnnotationFile", function(static, ...) {
  byChipType(static, ...);
}, static=TRUE, protected=TRUE)


setMethodS3("fromChipType", "AromaChipTypeAnnotationFile", function(static, ...) {
  byChipType(static, ...);
}, static=TRUE, protected=TRUE)





###########################################################################/**
# @RdocMethod findByChipType
#
# @title "Locates an annotation file by its chip type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chipType}{A @character string.}
#  \item{tags}{An optional @character @vector of tags.}
#  \item{...}{Additional arguments.}
# }
#
# \value{
#  Returns the pathname (as a @character string) to the first annotation
#  chip type file found.  If no one was found, @NULL is returned.
# }
#
# @author
#
# \seealso{
#   @seemethod "byChipType".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("findByChipType", "AromaChipTypeAnnotationFile", abstract=TRUE);





###########################################################################/**
# @RdocMethod getHeader
#
# @title "Gets the header of the annotation file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @list structure.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getHeader", "AromaChipTypeAnnotationFile", abstract=TRUE);


setMethodS3("getPlatform", "AromaChipTypeAnnotationFile", abstract=TRUE);


setMethodS3("getChipType", "AromaChipTypeAnnotationFile", abstract=TRUE);



############################################################################
# HISTORY:
# 2008-05-09
# o Created from AffymetrixCdfFile.R.
############################################################################
