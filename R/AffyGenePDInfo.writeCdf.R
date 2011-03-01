###########################################################################/** 
# @set "class=AffyGenePDInfo"
# @RdocMethod writeCdf
#
# @title "Generates an Affymetrix CDF file from a Platform Design (PD) package"
#
# \description{
#   @get "title".
#   Platform Design (PD) packages are also known as "pdInfo" packages.
# }
#
# @synopsis
#
# \arguments{
#  \item{tags}{An optional @character @vector of tags to be added to the CDF
#    filename.}
#  \item{path}{The output path where the CDF file is written.
#    If @NULL (default), then it is written to the corresponding
#    \code{annotationData/chipTypes/<chipType>/} directory.}
#  \item{overwrite}{If @TRUE, an existing CDF is overwritten, otherwise
#    an exception is thrown.}
#  \item{verbose}{A @logical or @see "R.utils::Verbose".}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) the pathname to CDF written.
# }
#
# \details{
#   The formal chip type of the CDF is inferred from the AffyGenePDInfo package.
#   The filename of the CDF is constructed from the chip type and any optional
#   tags.
#   To minimize the risk for a corrupt CDF file, the creation of the file
#   is atomic by first writing to a temporary file which is then renamed.
# }
#
# \section{Limitations}{
#   The information available in the PD package is limited and does
#   not contain all information needed to populate a CDF file.
#   In order to workaround these limitations, certain CDF entries
#   are set to predefined/hardwired values.
#   The 'pbase' and 'tbase' entries of the generated CDF file is
#   hardwired to "T" and "A", respectively.  Likewise, the 'groupdirection'
#   entry is hardwired to "sense".
# }
#
# \author{
#   Henrik Bengtsson, adopted from \code{pdInfo2Cdf()} written by
#   Samuel Wuest and Mark Robinson.
# }
#
# @keyword internal
#*/########################################################################### 
setMethodS3("writeCdf", "AffyGenePDInfo", function(this, tags=NULL, path=NULL, overwrite=FALSE, verbose=TRUE, ...) {
  require("affxparser") || throw("Package not loaded: affxparser");
  require("pdInfoBuilder") || throw("Package not loaded: pdInfoBuilder");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pmFeature2List <- function(u) {
    nr <- nrow(u);
    o <- order(u$atom);
    v <- 0:(nr-1);
    id <- u$fsetid[1];
    g <- list(list(x=u$x[o], y=u$y[o], pbase=rep("T", nr), tbase=rep("A", nr),
              atom=v, indexpos=v, groupdirection="sense", natoms=nr,
              ncellsperatom=1));
    names(g) <- id;
    list(unittype=1, unitdirection=1, groups=g, natoms=nr, ncells=nr,
         ncellsperatom=1, unitnumber=id);
  } # pmFeature2List()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- paste(tags, collapse=",");
    tags <- gsub(",,", ",", tags);
  }

  # Argument 'path':
  if (!is.null(path)) {
    path <- Arguments$getWritablePath(path);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);


  # Platform Design object
  pd <- this;

  verbose && enter(verbose, "Generating CDF file from Platform Design (PD) package");
  pkgName <- pd@annotation;
  verbose && cat(verbose, "Platform Design (PD) package: ", pkgName);

  # Infer the chip type
  pkgInfo <- packageDescription(pkgName);
  title <- pkgInfo$Title;
  chipType <- gsub(".* ", "", title);
  rm(pkgInfo, title);

  # Chip type package name
  pkgNameT <- cleanPlatformName(chipType);

  # Sanity check
  stopifnot(pkgNameT == pkgName);
  rm(pkgNameT);

  
  if (is.null(path)) {
    path <- file.path("annotationData", "chipTypes", chipType);
    path <- Arguments$getWritablePath(path);
  }

  verbose && cat(verbose, "Output path: ", path);

  chipTypeF <- paste(c(chipType, tags), collapse=",");
  filename <- sprintf("%s.cdf", chipTypeF);
  verbose && cat(verbose, "Filename: ", filename);

  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=!overwrite);
  verbose && cat(verbose, "Pathname to generated CDF: ", pathname);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve chip type dimensions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dim <- geometry(pd);
  nrows <- dim[1];
  ncols <- dim[2];
  rm(dim);

  nbrOfCells <- nrows*ncols;
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Tags: ", tags);
  verbose && printf(verbose, "Chip type dimensions: %dx%d\n", nrows, ncols);
  verbose && cat(verbose, "Total number of cells (probes): ", nbrOfCells);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving information from PD package
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Querying Platform Design database");
  ff <- dbGetQuery(db(pd), "select * from pmfeature");
  rm(pd);  # Not needed anymore
  verbose && str(verbose, ff);
  nbrOfPdCells <- nrow(ff);
  verbose && printf(verbose, "Number of cells (probes) in PD database: %d (%.2f%%) of %d\n",
                    nbrOfPdCells, 100*nbrOfPdCells/nbrOfCells, nbrOfCells);
  verbose && exit(verbose);

  verbose && enter(verbose, "Creating list from query table");
  # three 3 lines speed up the splitting ...
  ffs <- split(ff, substr(ff$fsetid, 1,4));
  ffs <- lapply(ffs, FUN=function(u) split(u, u$fsetid));
  ffs <- unlist(ffs, recursive=FALSE);
  names(ffs) <- substr(names(ffs), 6, nchar(names(ffs)));
  nbrOfUnits <- length(ffs);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);
  verbose && printf(verbose, "Average number of cells per units: %.2f\n", nbrOfPdCells/nbrOfUnits);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up CDF tree structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up CDF tree structure");

  verbose && enter(verbose, "Setting CDF header");
  newCdfHeader <- list(ncols=ncols, nrows=nrows, nunits=nbrOfUnits, 
                       nqcunits=0, refseq="", chiptype=chipType, 
                       filename=pathname, rows=nrows, 
                       cols=ncols, probesets=nbrOfUnits, 
                       qcprobesets=0, reference="");
  verbose && exit(verbose);

  verbose && enter(verbose, "Setting up CDF units");
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);
  newCdfList <- lapply(ffs, FUN=pmFeature2List);
  rm(ffs);  # Not needed anymore
  verbose && exit(verbose);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Writing CDF to file (binary format)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Writing (binary) CDF file");
  pathname <- newCdfHeader$filename;
  verbose && cat(verbose, "Pathname: ", pathname);

  # Write to a temporary file
  pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

  res <- affxparser::writeCdf(pathnameT, cdfheader=newCdfHeader, cdf=newCdfList,
                  cdfqc=NULL, verbose=verbose, overwrite=overwrite);

  # Rename temporary file
  pathname <- popTemporaryFile(pathnameT, verbose=verbose);

  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(res);
}) # writeCdf()



############################################################################
# HISTORY:
# 2011-01-09 [HB]
# o Added writeCdf() for AffyGenePDInfo, which replaces pdInfo2Cdf().
#   An auxillary CEL file is no longer needed to create a CDF from
#   an PDInfo package.  Moreover, contrary pdInfo2Cdf(), the generated
#   CDF now gets a correct/formal Affymetrix chip type.
# 
# Below history is for pdInfo2Cdf():
#
# 2010-12-04 [HB]
# o Added more verbose output.
# o DOCUMENTATION: Added more Rd documentation.
# o BUG FIX: Local variable 'pdName' of pdInfo2Cdf() was used before it
#   was defined.  Thanks to Guido Hooiveld at the Wageningen University, 
#   Netherlands, for reporting this.
# 2010-05-20 [HB]
# o Renamed PdInfo2Cdf() to pdInfo2Cdf().  Keeping old one for backward
#   compatibility for a while.
# 2010-05-19 [HB]
# o BUG FIX: PdInfo2Cdf() would write dimension (rows,rows) in the CDF 
#   header instead of (rows,cols).  Thanks Kasper Daniel Hansen for
#   reporting this.
# 2009-10-16 [HB]
# o MEMORY OPTIMIZATION: Cleaning out non needed objects sooner.
# o Added verbose a'la R.utils.
# o Added some validation of arguments.
# o Tidied up the code structure.
# 2009-01-13 [MR]
# o Added. "This script has been written to generate a .cdf-file from an 
#   "pd.XXXX" package, such as those build with pdInfoBuilder.
#   The original was written by Samuel Wuest, modified by Mark Robinson 
#   (around 12 Jan 2009) to be generic."
# 2008-??-?? [SW]
# o Created by Samuel Wuest.
############################################################################
