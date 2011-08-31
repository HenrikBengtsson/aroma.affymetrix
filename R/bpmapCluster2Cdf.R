###########################################################################/**
# @RdocDefault bpmapCluster2Cdf
#
# @title "Creates a CDF from tiling-array BPMAP file"
#
# \description{
#   @get "title".\cr
#
#   \emph{
#    NOTE: This method applies only to Affymetrix tiling arrays!
#    Furthermore, it is likely to be more useful for promoter tiling arrays
#    and less so for whole-genome tiling arrays.
#   }
# }
#
# @synopsis
#
# \arguments{
#  \item{pathname}{The pathname to the BPMAP file.}
#  \item{chipType, tags}{The chip type and optional tags of the CDF to
#    be written.}
#  \item{maxProbeDistance}{A positive @integer specifying the maximum
#    genomic distance (in basepairs) allowed between two probes in order
#    to "cluster" those two probes into the same CDF units.  Whenever the
#    distance is greater, the two probes end up in two different CDF units.}
#  \item{minNbrOfProbes}{A positive @integer specifying the minimum number
#    of probes required in a CDF unit.  If fewer, those probes are dropped.}
#  \item{...}{Not used.}
#  \item{rows, cols}{Two (optional) positive @integers.
#     If @NULL, optimal values are inferred auotmatically.}
#  \item{groupName}{A @character string.}
#  \item{field}{A @character string.}
#  \item{stringRemove}{An (optional) regular expression.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns (invisibly) a the pathname of the created CDF file.
#  The created CDF is written to the current directory.
# }
#
# \details{
#   This method applies only to Affymetrix tiling arrays.  It is likely
#   to be useful for promoter tiling arrays and less so for whole-genome
#   tiling arrays.
# }
#
# \author{
#   Henrik Bengtsson adopted from Mark Robinson standalone/online version
#   as of July 11, 2011.
# }
# 
# @keyword "internal"
#*/###########################################################################
setMethodS3("bpmapCluster2Cdf", "default", function(pathname, chipType, tags=NULL, maxProbeDistance=3000, minNbrOfProbes=30, ..., rows=NULL, cols=NULL, groupName=gsub("_.*", "", chipType), field="fullname", stringRemove=sprintf("%s:.*;", groupName), verbose=-10) {
  require("affxparser") || throw("Package not loaded: affxparser");
  require("R.utils") || throw("Package not loaded: R.utils");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  bpmapUnit2df <- function(u) {
    o <- order(u[["startpos"]]);
    mmx <- u[["mmx"]];
    mmy <- u[["mmy"]];
    mmx <- if (all(mmx == 0L) | is.null(mmx)) 0L else mmx;
    mmy <- if (all(mmy == 0L) | is.null(mmy)) 0L else mmy;

    seqInfo <- u$seqInfo;
    res <- data.frame(seqname=seqInfo[[field]], groupname=seqInfo$groupname, u[c("pmx", "pmy")], mmx=mmx, mmy=mmy, u[c("probeseq", "strand", "startpos", "matchscore")], stringsAsFactors=FALSE);
    res <- res[o,,drop=FALSE];
    res;
  } # bpmapUnit2df()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validating arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'filename':
  pathname <- Arguments$getReadablePathname(pathname);

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
  }

  # Argument 'minNbrOfProbes':
  minNbrOfProbes <- Arguments$getInteger(minNbrOfProbes, range=c(1,Inf));

  # Argument 'maxProbeDistance':
  maxProbeDistance <- Arguments$getInteger(maxProbeDistance, range=c(1,Inf));

  # Argument 'rows':
  if (!is.null(rows)) {
    rows <- Arguments$getInteger(rows, range=c(1,Inf));
  }

  # Argument 'cols':
  if (!is.null(cols)) {
    cols <- Arguments$getInteger(cols, range=c(1,Inf));
  }

  # Argument 'groupName':
  groupName <- Arguments$getCharacter(groupName);

  # Argument 'field':
  field <- Arguments$getCharacter(field);

  # Argument 'stringRemove':
  if (!is.null(stringRemove)) {
    stringRemove <- Arguments$getCharacter(stringRemove);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);



  verbose && enter(verbose, "Generating CDF from BPMAP");

  path <- ".";  
  path <- Arguments$getWritablePath(path);
  verbose && cat(verbose, "Directory where CDF will be written: ", path);

  fullname <- paste(c(chipType, tags), collapse=",");
  verbose && cat(verbose, "Full name: ", fullname);

  cdfFilename <- sprintf("%s.cdf", fullname);
  cdfFilename <- Arguments$getFilename(cdfFilename);
  verbose && cat(verbose, "CDF pathname: ", cdfFilename);
  cdfPathname <- Arguments$getWritablePathname(cdfFilename, mustNotExist=TRUE);

  ppsPathname <- sprintf("%s.pps", fullname);
  ppsPathname <- Arguments$getFilename(ppsPathname);
  verbose && cat(verbose, "PPS pathname: ", ppsPathname);
  ppsPathname <- Arguments$getWritablePathname(ppsPathname, mustNotExist=TRUE);

  verbose && cat(verbose, "Argument 'groupName': ", groupName);
  verbose && cat(verbose, "Argument 'stringRemove': ", stringRemove);

  verbose && enter(verbose, "Reading BPMAP file");
  verbose && cat(verbose, "Source pathname: ", pathname);
  bpmapList <- readBpmap(pathname, readMatchScore=TRUE);
  verbose && cat(verbose, "Number of BPMAP units: ", length(bpmapList));
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting X/Y locations");
  bpmapdfList <- lapply(bpmapList, FUN=bpmapUnit2df);
  verbose && exit(verbose);

  rm(bpmapList); # Not needed anymore


  # Infer number of CDF rows and columns, if missing.
  if (is.null(rows)) {
    rows <- max(sapply(bpmapdfList, FUN=function(u) max(c(u$mmx,u$pmx))));
    verbose && printf(verbose, "NB: 'rows' of CDF are being set as %d. If this is not correct, stop now and specify 'rows' argument.\n", rows);
  }

  if (is.null(cols)) {
    cols <- max(sapply(bpmapdfList, FUN=function(u) max(c(u$mmy,u$pmy))));
    verbose && printf(verbose, "NB: 'cols' of CDF are being set as %d. If this is not correct, stop now and specify 'col' argument.\n", cols);
  }

  verbose && enter(verbose, "Counting the number of CDF units");
  nbrOfUnits <- 0L;
  for (ii in seq(along=bpmapdfList)) {
    name <- names(bpmapdfList)[ii];
    verbose && enter(verbose, sprintf("Item #%d ('%s') of %d", ii, name, length(bpmapdfList)));
    bpmapdf <- bpmapdfList[[ii]];
    np <- nrow(bpmapdf);

    sp <- bpmapdf$startpos;
    if (all(sp > 0L) & bpmapdf$groupname[1] == groupName) {
      d <- diff(sp);
      w <- whichVector(d > maxProbeDistance);
      ends <- c(w, np);
      starts <- c(1L, w+1L);
      k <- whichVector((ends-starts) > minNbrOfProbes);
      verbose && cat(verbose, length(k), " ROIs for ", name, ".");

      nbrOfUnits <- nbrOfUnits + length(k);
    } else {
      # keep all probes
      verbose && cat(verbose, "Skipping all ", np, " probes.");
    }
    verbose && exit(verbose);
  } # for (ii ...)
  verbose && cat(verbose, "Number of CDF units: ", nbrOfUnits);

  # Sanity check
  stopifnot(nbrOfUnits > 0);
  verbose && exit(verbose);

    
  verbose && enter(verbose, "Creating CDF tree structure for ", length(bpmapdfList), " BPMAP units");

  # Allocate the CDF tree structure
  cdfList <- vector("list", length=nbrOfUnits);
  naValue <- as.character(NA);
  unitNames <- rep(naValue, times=length(cdfList));

  startps <- vector("list", length=nbrOfUnits);
  e <- vector("list", length=1L);
  uu <- 1L;
  for (ii in seq(along=bpmapdfList)) {
    name <- names(bpmapdfList)[ii];
    verbose && enter(verbose, sprintf("Item #%d ('%s') of %d", ii, name, length(bpmapdfList)));

    # Access ones
    bpmapdf <- bpmapdfList[[ii]];
    np <- nrow(bpmapdf);

    ch <- bpmapdf$seqname[1];
    if (!is.null(stringRemove)) {
      ch <- gsub(stringRemove, "", ch);
    }

    #if (all(bpmapdf$mmx == 0L) && all(bpmapdf$startpos > 0L)) {
    sp <- bpmapdf$startpos;
    if (all(sp > 0L) & bpmapdf$groupname[1] == groupName) {
      d <- diff(sp);
      w <- whichVector(d > maxProbeDistance);
      ends <- c(w, np);
      starts <- c(1L, w+1L);
      k <- whichVector((ends-starts) > minNbrOfProbes);
      verbose && cat(verbose, length(k), " ROIs for ", name, ".");

      # Access ones
      pmx <- bpmapdf$pmx;
      pmy <- bpmapdf$pmy;
      mmx <- bpmapdf$mmx;
      mmy <- bpmapdf$mmy;

      for (jj in seq(along=k)) {
        w <- starts[k[jj]]:ends[k[jj]];
        np <- length(w);

        atom <- 0L:(np-1L);
        indexpos <- 0L:(np-1L);

        if (all(mmx == 0)) {
          # PM only
          e[[1]] <- list(x=pmx[w], y=pmy[w], pbase=rep("A", times=np), tbase=rep("T", times=np), atom=atom, indexpos=indexpos, groupdirection="sense", natoms=np, ncellsperatom=1L); 
        } else {
          # PM+MM
          e[[1]] <- list(x=c(pmx[w],mmx[w]), y=c(pmy[w],mmy[w]), pbase=rep("A", times=2*np), tbase=rep(c("T","A"), each=np), atom=rep(atom, times=2), indexpos=rep(indexpos, times=2), groupdirection="sense", natoms=np, ncellsperatom=2L); 
        }

        names(e) <- paste(ch, "FROM", sp[starts[k[jj]]], "TO", sp[ends[k[jj]]], sep="");
        na <- sum(unlist(sapply(e, FUN=function(u) u$natoms), use.names=FALSE));
        nc <- sum(unlist(sapply(e, FUN=function(u) u$natoms*u$ncellsperatom), use.names=FALSE));

        cdfList[[uu]] <- list(unittype=1L, unitdirection=1L, groups=e, natoms=na, ncells=nc, ncellsperatom=nc/na, unitnumber=ii);
        startps[[uu]] <- sp[w];
        unitNames[uu] <- names(e);

        # Next CDF unit
        uu <- uu + 1L;
      } # for (jj ...)
    } else {
      # keep all probes
      verbose && cat(verbose, "Skipping.");
    }

    verbose && exit(verbose);
  } # for (ii ...)
  verbose && exit(verbose);

  names(cdfList) <- unitNames;

  # Sanity check
  stopifnot(length(cdfList) == nbrOfUnits);


  verbose && enter(verbose, "Writing PPS file");
  verbose && cat(verbose, "Output pathname: ", ppsPathname);
  saveObject(startps, file=ppsPathname);
  rm(startps); # Not needed anymore
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing CDF file");
  verbose && cat(verbose, "Output pathname: ", cdfPathname);

  cdfHeader <- list(probesets=nbrOfUnits, qcprobesets=0, reference="", chiptype=chipType, filename=cdfPathname, nqcunits=0, nunits=nbrOfUnits, rows=rows, cols=cols, refseq="", nrows=rows, ncols=cols);
  verbose && cat(verbose, "CDF header:");
  verbose && str(verbose, cdfHeader);

  # Write to a temporary file
  pathnameT <- pushTemporaryFile(cdfPathname, verbose=verbose);
 
  writeCdf(pathnameT, cdfheader=cdfHeader, cdf=cdfList, cdfqc=NULL, overwrite=TRUE, verbose=verbose);

  # Rename temporary file
  pathname <- popTemporaryFile(pathnameT, verbose=verbose);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(pathname);
}) # bpmapCluster2Cdf()


############################################################################
# HISTORY:
# 2011-08-31 [HB]
# o All arguments after '...' except 'verbose' may be dropped in a future
#   release, especially 'rows' and 'cols'.
# o CLARIFICATION: Renamed argument 'nProbes' to 'minNbrOfProbes'.
# 2011-08-30 [HB]
# o Now bpmapCluster2Cdf() returns the pathname to the created CDF.
# 2011-08-29 [HB]
# o Replaced argument 'cdfName' with 'chipType' and 'tags'.  This 
#   simplifies adding custom tags to the CDF.
#   It will also allow us (later) to write the CDF to the correct
#   annotationData/ directory.
# o ROBUSTNESS: Now the writing of the CDF file is atomic by first writing
#   to a temporary file which is then renamed.
# o GENERALIZATION/ROBUSTNESS: Now the method precount the number of
#   CDF units so it can allocate a CDF tree structure of the correct size.
# o GENERALIZATION: Now the default value of argument 'groupName' is 
#   more general.  Before it was hardcoded to "Hs".
# o GENERALIZATION: Now the default value of argument 'stringRemove' is 
#   more general.  Before it was hardcoded to a particular BPMAP file.
# o DOCUMENTATION: Added some basic Rdoc documentation.
# o Added more verbose progress output.
# o CLEANUP: Allocated variable 'll' was never used.
# 2011-07-11 [HB]
# o CLEANUP: Code and argument cleanup.
# o ROBUSTNESS: Using '&&' instead of '&' if-statement.
# o ROBUSTNESS: Added more argument validation.
# o ROBUSTNESS: Added protection against overwriting existing CDFs.
# o ROBUSTNESS: Added tests for valid CDF filenames.
# 2008-11-28 [HB]
# o CLEANUP: Tidying up code.
# o Now using R.utils::Verbose statements.
# 2008-11-xx [MR]
# o Created.
############################################################################
