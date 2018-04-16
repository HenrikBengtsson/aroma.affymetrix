# @author "HB"
setConstructorS3("AffymetrixTabularFile", function(...) {
  extend(TabularTextFile(...), c("AffymetrixTabularFile",
              uses("AromaPlatformInterface"), uses("FileCacheKeyInterface"))
  )
})


setMethodS3("translateColumnNames", "AffymetrixTabularFile", function(this, names, ...) {
  # Convert 'FOO_BAR_DOO' and 'FOO.BAR.DOO' to 'foo bar doo'?
  if (any(regexpr("[_.]", names) != -1)) {
    names <- tolower(gsub("[_.]", " ", names))
  }

  # Finally, convert 'Foo bar Doo' to 'fooBarDoo'
  names <- toCamelCase(names)

  names
}, protected=TRUE)


setMethodS3("findByChipType", "AffymetrixTabularFile", function(static, chipType, tags=NULL, pattern=NULL, ...) {
  if (is.null(pattern)) {
    name <- paste(c(chipType, tags), collapse=",")
    pattern <- sprintf("^%s.*[.]...$", name)
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- findAnnotationDataByChipType(chipType, pattern, ...)
  pathname
}, static=TRUE, protected=TRUE)



setMethodS3("byChipType", "AffymetrixTabularFile", function(static, chipType, tags=NULL, ...) {
  # Search for the genome information file
  pathname <- findByChipType(static, chipType, tags=tags, ...)
  if (is.null(pathname))
    throw("Failed to located Affymetrix tabular file: ", chipType)
  newInstance(static, pathname, ...)
}, static=TRUE)
