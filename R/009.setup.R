.setupAromaAffymetrix <- function(pkg, ...) {
  # To please R CMD check
  ns <- getNamespace("aroma.core")
  .requireBiocPackage <- get(".requireBiocPackage", mode="function", envir=ns)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # None at the moment.

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Bioconductor packages aroma.light and affxparser
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # require("aroma.light") - install if missing
  .requireBiocPackage("aroma.light", neededBy=getName(pkg))

  # require("affxparser") - install if missing
  .requireBiocPackage("affxparser", neededBy=getName(pkg))

  # Make sure 'affxparser' is after 'aroma.affymetrix' on the search path
  from <- "package:affxparser"
  to <- "package:aroma.affymetrix"
  fromIdx <- match(from, search())
  toIdx <- match(to, search())
  if (all(is.finite(c(fromIdx, toIdx))) && fromIdx < toIdx) {
    moveInSearchPath(from=from, to=to, where="after")
  }

  # Add custom findCdf() function to affxparser.  This is need to be
  # able to locate CDFs in annotationData/chipTypes/<chipType>/.
  setCustomFindCdf(function(...) {
    AffymetrixCdfFile$findByChipType(..., .useAffxparser=FALSE)
  })


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Package settings (settings might change)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # This code will update the settings according to the default ones.
  updateSettings(pkg)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fix the search path every time a package is loaded
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setHook("base::library:onLoad", function(...) {
    # Fix the search path
    pkgs <- fixSearchPath(aroma.affymetrix)
    if (length(pkgs) > 0) {
      warning("Packages reordered in search path: ",
                                            paste(pkgs, collapse=", "))
    }
  }, action="append")
} # .setupAromaAffymetrix()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# affxparser related
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Instead of asking users to write affxparser::writeCdf() when
# aroma.affymetrix is loaded...
setMethodS3("writeCdf", "default", function(...) {
  ns <- loadNamespace("affxparser")
  `affxparser::writeCdf` <- get("writeCdf", envir=ns, mode="function")
  `affxparser::writeCdf`(...)
})
