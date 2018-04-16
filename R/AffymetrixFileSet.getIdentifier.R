setMethodS3("getIdentifier", "AffymetrixFileSet", function(this, ...) {
  path <- getPath(this)
  res <- NULL
  for (kk in 1:3) {
    pathname <- file.path(path, "IDENTIFIER")
    if (isFile(pathname)) {
      res <- readLines(pathname)
      # Remove comments
      res <- trim(gsub("#.*", "", trim(res)))
      # Remove empty lines
      res <- res[nzchar(res)]
      break
    }
    path <- dirname(path)
  }

  if (!is.null(res)) {
    res <- getChecksum(list(res))
  }

  res
}, private=TRUE)
