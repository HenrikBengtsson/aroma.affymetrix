###########################################################################/**
# @set "class=AffymetrixCelFile"
# @RdocMethod createFrom
#
# @title "Creates a CEL file using another as a template"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filename, path}{The filename and path of to the CEL 
#     file to be created.}
#  \item{version}{The file-format version of the CEL file to be created.}
#  \item{methods}{If \code{"copy"}, the new file is created as a copy of the
#     template file.  If \code{"create"}, the new file is created from
#     scratch from the template file.}
#  \item{clear}{If @TRUE, the fields of the CEL file are cleared (zeroed),
#     otherwise they contain the same information as the source file.}
#  \item{...}{Not used.}
#  \item{verbose}{See "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "AffymetrixCelFile" reference to the new CEL file.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("createFrom", "AffymetrixCelFile", function(this, filename, path=NULL, overwrite=FALSE, skip=!overwrite, version=c("4", "3"), methods=c("copy", "create"), clear=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Rename lower-case *.cel to *.CEL, if that is the case.  Old versions
  # of the package generated lower-case CEL files. /HB 2007-08-09
  if (regexpr("[.]cel$", pathname) != -1) {
    pathname <- AffymetrixFile$renameToUpperCaseExt(pathname); 
  }

  pathname <- Arguments$getWritablePathname(pathname, 
                                          mustNotExist=(!overwrite && !skip));

  # Argument 'version':
  version <- match.arg(version);

  # Argument 'methods':
  if (!all(methods %in% c("copy", "create"))) {
    throw("Unknown value of argument 'methods': ", 
                                              paste(methods, collapse=", "));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Get the CDF of the template CEL file
  cdf <- getCdf(this);

  verbose && enter(verbose, "Creating CEL file");
  verbose && cat(verbose, "Chip type: ", getChipType(cdf));
  verbose && cat(verbose, "Pathname: ", pathname);

  # Nothing to do?
  if (skip && isFile(pathname)) {
    verbose && cat(verbose, "Returning already existing file.");
    res <- newInstance(this, pathname);
    ver <- getHeader(res)$version;
    if (ver != version) {
      throw("Cannot not retrieve CEL file of version ", version, 
         ". The existing CEL file has version ", ver, ": ", pathname);
    }
    setCdf(res, cdf);
    verbose && exit(verbose);
    return(res);
  }


  # First create/copy to a temporary file, then rename
  tmpPathname <- paste(pathname, "tmp", sep=".");
  if (isFile(tmpPathname)) {
    throw("Cannot create CEL file. Temporary file already exists: ", tmpPathname);
  }


  msgs <- list();
  res <- NULL;
  for (method in methods) {
    verbose && enter(verbose, "Method '", method, "'");

    if (method == "copy") {
      # Check version of template CEL file
      ver <- getHeader(this)$version;
      if (ver != version) {
        msgs[[method]] <- paste("Cannot create CEL file of version ", version,
           " (", pathname, "). Template CEL file is of version ", 
           ver, ": ", getPathname(this), sep="");
        verbose && cat(verbose, msgs[[method]]);
        verbose && exit(verbose);
        next;
      }

      # 1. Create a temporary file    
      res <- copyTo(this, filename=tmpPathname, path=NULL, 
                                             verbose=less(verbose));
#      verbose && cat(verbose, "Temporary file:");
#      verbose && print(verbose, res);
    
      # 2. Update the temporary file
      if (clear) {
        clearData(res, ..., .forSure=TRUE, verbose=less(verbose));
      }

      # 3. Rename the temporary file
      renameTo(res, filename=pathname, verbose=less(verbose));

      # Break out of the methods loop
      verbose && exit(verbose);
      break;
    }
  
    if (method == "create") {
      if (version != "4") {
        msgs[[method]] <- paste(
          "Can only create binary CEL files (version 4) from scratch, ",
          "not files of version ", version, ": ", pathname, sep="");
        verbose && cat(verbose, msgs[[method]]);
        verbose && exit(verbose);
        next;
      }
  
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Setting up CEL header
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      fullname <- getFullName(this);
      celHeader <- cdfHeaderToCelHeader(getHeader(cdf), sampleName=fullname);
    
      # Add some extra information about what the CEL file is for
      params <- c(Descripion=sprintf("This CEL file contains data saved by the aroma.affymetrix v%s package.", getVersion(aroma.affymetrix)));
      parameters <- gsub(" ", "_", params);
      names(parameters) <- names(params);
      parameters <- paste(names(parameters), parameters, sep=":");
      parameters <- paste(parameters, collapse=";");
      parameters <- paste(celHeader$parameters, parameters, "", sep=";");
      parameters <- gsub(";;", ";", parameters);
      parameters <- gsub(";$", "", parameters);
      celHeader$parameters <- parameters;
      rm(params, parameters);
    
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Creating empty CEL file
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # 1. Create a temporary file    
      verbose && enter(verbose, "Creating an empty temporary CEL file");
      createCel(tmpPathname, header=celHeader, overwrite=overwrite, ...,
                                                        verbose=less(verbose));
      rm(celHeader);
      verbose && exit(verbose);

      # 2. Update the temporary file
      if (!clear) {
        verbose && enter(verbose, "Copying CEL data");
        cells <- seq(length=nbrOfCells(this));
        lapplyInChunks(cells, function(cells) {
          verbose && enter(verbose, "Reading subset of data from source CEL file");
          data <- readCel(getPathname(this), indices=cells, readIntensities=TRUE, readStdvs=TRUE, readPixels=TRUE);
          verbose && str(verbose, data, level=-50);
          verbose && printf(verbose, "RAM: %.2fMB\n", object.size(data)/1024^2, level=-40);
          verbose && exit(verbose);
          gc <- gc();

          verbose && enter(verbose, "Writing data to new CEL file");
          updateCel(tmpPathname, indices=cells, intensities=data);
          verbose && exit(verbose);

          rm(data);
          gc <- gc();
          verbose && print(verbose, gc);
        }, chunkSize=1e6, verbose=verbose);
        rm(cells);

        verbose && exit(verbose);
      }
  
      # 3. Rename the temporary file
      verbose && enter(verbose, "Renaming file");
      res <- file.rename(tmpPathname, pathname);
      if (!res) {
        throw("Failed to rename temporary file: ", 
                                   tmpPathname, " -> ", pathname);
      }
      verbose && exit(verbose);
      
      res <- newInstance(this, pathname);

      # Break out of the methods loop
      verbose && exit(verbose);
      break;
    }

    verbose && exit(verbose);
  } # for (method ...)

  if (is.null(res)) {
    msgs <- unlist(msgs);
    msgs <- paste(names(msgs), msgs, sep=": ");
    msgs <- paste(msgs, collapse="; ");
    msg <- paste("Failed to create CEL file. Error messages: ", msgs);
    verbose && cat(verbose, "Error: ", msg);
    throw(msg);
  }

  # Make sure the CDF is carried down
  setCdf(res, cdf);

  verbose && print(verbose, res);

  verbose && exit(verbose);

  res;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2008-12-12
# o createFrom() for AffymetrixCelFile had verbose=TRUE as default.
# 2007-09-16
# o Made createFrom() fail safe, that is, the destination file will not
#   be created unless it is complete.  This is done by first create/writing
#   to a temporary file which is then renamed at the end.
# o Now createFrom() is using copyTo() instead of obsolete copyFile().
# 2007-08-16
# o Renamed argument 'method' in createFrom() of AffymetrixCelSet to
#   'methods' to reflect the fact that it can take multiple values.
#   The default is not that createFrom() first tries to copy, and if that
#   fails/is not possible, it tries to create the CEL file from scratch.
#   For instance, calling createFrom() on an binary XDA CEL file, will copy
#   it, but calling it on a Calvin CEL file, will create a binary XDA CEL
#   file (because aroma.affymetrix/affxparser do not support updating 
#   Calvin CEL file so we do not support "creating" them.
# 2007-03-24
# o Argument 'method' still defaults to "copy" altough "create" is much
#   faster.  Argument 'clear' now defaults to FALSE.  The reason for this 
#   is that we want in most cases want to "inherit" some of the fields
#   from the source, e.g. we don't want to change stdvs in normalization.
# o Added code to create CEL files from scratch by adopting the code in
#   createBlankSet() of AffymetrixCelSet.
# o Moved to its own file.
# 2007-03-23
# o Added argument 'version' to createFrom(), which defaults to CEL file
#   format version 4.  Now it is asserted that the template CEL file is
#   of the same version.
# 2006-08-27
# o Added createFrom() which utilizes new functions copyFile() and 
#   clearData(). It is also no longer static. This is more generic and 
#   cleaner.  The new clearData() does also not require the CDF file 
#   (in case that should be missing).
############################################################################
