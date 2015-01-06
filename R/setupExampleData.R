###########################################################################/**
# @set class=AromaAffymetrix
# @RdocMethod setupExampleData
#
# @title "Setups example data in the current directory"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dirs}{A @character @vector specifying which directories to setup.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) @TRUE if all requested data was installed,
#   otherwise @FALSE.
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("setupExampleData", "AromaAffymetrix", function(pkg, dirs=c("annotationData", "rawData"), mustWork=TRUE, validate=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dirs':
  dirs <- Arguments$getCharacters(dirs);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check needed packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (pkg in c("AffymetrixDataTestFiles", "affxparser")) {
    if (!isPackageInstalled(pkg)) {
      if (!mustWork) return(invisible(FALSE));
      throw("Package not installed: ", pkg);
    }
  }

  path0 <- system.file(package="AffymetrixDataTestFiles")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # annotationData/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathR <- "annotationData"
  if (is.element(pathR, dirs)) {
    chipType <- "HG-Focus"
    pathD <- file.path(pathR, "chipTypes", chipType)
    filename <- "HG-Focus.CDF";
    pathnameD <- Arguments$getReadablePathname(filename, path=pathD, mustExist=FALSE)

    # Missing?
    if (!isFile(pathnameD)) {
      pathnameD <- Arguments$getWritablePathname(pathnameD)
      pathS <- file.path(path0, pathD, "3.ASCII")
      pathnameS <- Arguments$getReadablePathname(filename, path=pathS)

      cdfS <- AffymetrixCdfFile(pathnameS)
      formatS <- getFileFormat(cdfS)
      if (regexpr("v4", formatS) != -1L) {
        # Link
        createLink(link=pathnameD, target=pathnameS)
      } else {
        # Convert source ASCII CDF to binary CDF
        .convertCdf(pathnameS, pathnameD, .validate=FALSE, verbose=TRUE)
      }

      # Sanity check
      if (validate) {
        cdfD <- AffymetrixCdfFile(pathnameD)
        formatD <- getFileFormat(cdfD)
        stopifnot(regexpr("v4", formatD) != -1L)
      }
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # rawData/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathR <- "rawData"
  if (is.element(pathR, dirs)) {
    dataset <- "FusionSDK_HG-Focus"
    chipType <- "HG-Focus"
    pathD <- file.path(pathR, dataset, chipType)
    pathD <- Arguments$getWritablePath(pathD)

    pathS <- file.path(path0, pathD)
    pathS <- file.path(pathS, "2.Calvin") # NOTE!

    pathnamesS <- list.files(path=pathS, pattern="[.]CEL$", full.names=TRUE)
    for (ii in seq_along(pathnamesS)) {
      pathnameS <- pathnamesS[ii]
      pathnameD <- file.path(pathD, basename(pathnameS))
      createLink(link=pathnameD, target=pathnameS, skip=TRUE)
    }

    # Sanity check
    if (validate && isDirectory("annotationData")) {
      cels <- AffymetrixCelSet$byName(dataset, chipType=chipType)
    }
  }

  invisible(TRUE);
}, protected=TRUE) # setupExampleData()


############################################################################
# HISTORY:
# 2014-06-24
# o Created from aroma.seq ditto.
############################################################################
