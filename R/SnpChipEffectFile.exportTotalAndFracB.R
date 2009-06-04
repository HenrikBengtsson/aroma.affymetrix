setMethodS3("exportTotalAndFracB", "SnpChipEffectFile", function(this, fields=c("total", "fracB"), fullname=gsub(",chipEffects", "", getFullName(this)), dataSet=NULL, path=NULL, rootPath="totalAndFracBData", ..., overwrite=FALSE, drop=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'fullname':
  fullname <- Arguments$getCharacter(fullname, length=c(1,1));

  # Argument 'field':
  fields <- match.arg(fields, several.ok=TRUE);

  # Arguments 'dataSet':
  if (is.null(dataSet)) {
    dataSet <- basename(getParent(getPath(this)));
  }
  dataSet <- Arguments$getCharacter(dataSet, length=c(1,1));

  # Arguments 'path':
  if (is.null(path)) {
    chipType <- getChipType(this, fullname=FALSE);
    path <- filePath(rootPath, dataSet, chipType);
    rm(chipType);
  }
  path <- Arguments$getWritablePath(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Exporting (total, fracB) data");
  verbose && cat(verbose, "Data set: ", dataSet);
  verbose && cat(verbose, "Fullname: ", fullname);
  verbose && cat(verbose, "Path: ", path);

  verbose && cat(verbose, "Fields: ", paste(fields, collapse=", "));

  cdf <- getCdf(this); 
  nbrOfUnits <- nbrOfUnits(cdf);
  platform <- getPlatform(this);
  chipType <- getChipType(this);
  chipType <- gsub(",monocell", "", chipType, fixed=TRUE);

  footer <- list(
    srcFile=list(
      srcDataSet=dataSet,
      srcChipType=getChipType(this),
      srcFullName=getFullName(this),
      srcChecksum=getChecksum(this)
    )
  );

  data <- NULL;

  asbList <- list();
  for (field in fields) {
    # Identify output class
    if (field == "total") {
      signalClass <- AromaUnitTotalCnBinaryFile;
    } else if (is.element(field, c("fracB", "freqB"))) {
      signalClass <- AromaUnitFracBCnBinaryFile;
    }
  
    verbose && enter(verbose, "Exporting ", class(this)[1], " as an ", getName(signalClass));
    verbose && cat(verbose, "Signal: ", field);
  
    # Generate output filename
    filename <- sprintf("%s,%s.asb", fullname, field);
    pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=FALSE); 
    verbose && cat(verbose, "Output pathname: ", pathname);
    if (isFile(pathname)) {
      if (!overwrite) {
        verbose && cat(verbose, "Output file already exists. Return that instead.");
        asb <- signalClass$fromFile(pathname);
        asbList[[field]] <- asb;
        verbose && exit(verbose);
        next;
      }
    }
  
    verbose && cat(verbose, "File footer:");
    verbose && str(verbose, footer);
  
    # Reading data
    if (is.null(data)) {
      verbose && enter(verbose, "Reading data");
      data <- extractTotalAndFracB(this, verbose=less(verbose, 5));
      # Backward compatibility
      colnames(data) <- gsub("freqB", "fracB", colnames(data), fixed=TRUE);
      verbose && str(verbose, data);
      verbose && exit(verbose);
    }

    verbose && cat(verbose, "Values:");
    values <- data[,field, drop=TRUE];
    verbose && str(verbose, values);

    verbose && enter(verbose, "Allocating output file");
    asb <- signalClass$allocate(filename=pathname, path=NULL, nbrOfRows=nbrOfUnits, platform=platform, chipType=chipType, footer=footer, overwrite=overwrite, verbose=less(verbose, 25));
    verbose && print(verbose, asb);
    verbose && exit(verbose);
  
    verbose && enter(verbose, "Writing data");
    asb[,1] <- values;
    verbose && exit(verbose);
    
    verbose && exit(verbose);

    asbList[[field]] <- asb;
  } # for (field ...)
  names(asbList) <- fields;
  rm(data);

  if (drop && length(asbList) == 1) {
    asbList <- asbList[[1]];
  }

  verbose && exit(verbose);

  invisible(asbList);
}, protected=TRUE) # exportTotalAndFracB()


setMethodS3("exportTotalAndFracB", "SnpChipEffectSet", function(this, fields=c("total", "fracB"), rootPath="totalAndFracBData", ..., drop=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fields':
  fields <- match.arg(fields, several.ok=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  signalClassList <- lapply(fields, FUN=function(field) {
    if (field == "total") {
      signalClass <- AromaUnitTotalCnBinarySet;
    } else if (is.element(field, c("fracB", "freqB"))) {
      signalClass <- AromaUnitFracBCnBinarySet;
    }
    signalClass;
  });

  names <- paste(sapply(signalClassList, FUN=getName), collapse=" and ");
  verbose && enter(verbose, "Exporting ", class(this)[1], " as ", names);

  dataSetName <- getFullName(this);
  chipType <- NULL;
  for (kk in seq(this)) {
    cf <- getFile(this, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, getName(cf), length(this)));
    asbList <- exportTotalAndFracB(cf, fields=fields, dataSet=dataSetName, ..., drop=FALSE, verbose=less(verbose, 1));
    if (is.null(chipType)) {
      chipType <- getChipType(asbList[[1]], fullname=FALSE);
    }
    verbose && print(verbose, asbList);
    rm(asbList);
    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);
  

  assList <- lapply(signalClassList, function(signalClass) {
    verbose && enter(verbose, "Setting up the ", getName(signalClass));
    ass <- NULL;
    tryCatch({
      ass <- signalClass$byName(dataSetName, chipType=chipType, paths=rootPath);
      verbose && print(verbose, ass);
    }, error = function(ex) {
    })
    verbose && exit(verbose);
    ass;
  });
  names(assList) <- fields;

  assList <- assList[!sapply(assList, is.null)];

  if (drop && length(assList) == 1) {
    assList <- assList[[1]];
  }

  invisible(assList);
}, protected=TRUE) # exportTotalAndFracB()


setMethodS3("exportTotalAndFracB", "CnChipEffectFile", function(this, fields=c("total", "fracB"), ...) {
  # Don't export fracB signals, if they are not available
  if (this$combineAlleles) {
    fields <- setdiff(fields, "fracB");
  }

  NextMethod("exportTotalAndFracB", this, fields=fields, ...);
})


setMethodS3("exportTotalAndFracB", "CnChipEffectSet", function(this, fields=c("total", "fracB"), ...) {
  # Don't export fracB signals, if they are not available
  if (getCombineAlleles(this)) {
    fields <- setdiff(fields, "fracB");
  }

  NextMethod("exportTotalAndFracB", this, fields=fields, ...);
})



setMethodS3("getAromaUnitTotalCnBinarySet", "default", function(this, ...) {
  exportTotalAndFracB(this, fields="total", ...);
})

setMethodS3("getAromaUnitFracBCnBinarySet", "default", function(this, ...) {
  exportTotalAndFracB(this, fields="fracB", ...);
})



############################################################################
# HISTORY:
# 2009-02-24
# o BUG FIX: exportTotalAndFracB() of SnpChipEffectFile return an empty
#   list for chip types with tags.
# 2009-02-22
# o Now exportTotalAndFracB() of CnChipEffect{File|Set} does not export
#   fracB signals if allele-specific chip effects do not exist.
# o exportTotalAndFracB() of SnpChipEffectFile would write the short
#   chip type in the file footer, not the full one.  This could lead to 
#   using the wrong annotation files etc.
# 2009-02-11
# o Now exported chip effect files no longer contains tag 'chipEffects'.
# o Renamed all methods.
# o Added argument 'rootPath'.
# 2008-09-10
# o Updated to be compatible with new aroma.core.
# 2008-07-30
# o Added getTotalAndFreqBSets() which is a more convenient name.
# 2008-06-25
# o Created.
############################################################################ 
