# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# BEGIN: AFFX specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setMethodS3("getListOfGenomeInformations", "ChromosomalModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Retrieving genome informations")
  cdfList <- getListOfCdfs(getSetTuple(this), ...)
  giList <- lapply(cdfList, FUN=getGenomeInformation, verbose=less(verbose))
  verbose && exit(verbose)

  giList
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# END: AFFX specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
