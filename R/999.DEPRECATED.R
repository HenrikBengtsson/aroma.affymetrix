# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# [() and [[() should be used to extract files (and nothing else)
#
# 2012-11-20
# o CLEANUP: Deprecated "[" and "[[" for AffymetrixCelFile,
#   AffymetrixCelSet CnagCfhFile, and CnagCfhSet, because they
#   should be used to subset files (and not units).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("[", "AffymetrixCelFile", function(this, units=NULL, drop=FALSE) {
  .Deprecated("readUnits");
  data <- readUnits(this, units=units);
  if (drop && length(data) == 1)
    data <- data[[1]];
  data;
}, protected=TRUE, deprecated=TRUE)

setMethodS3("[[", "AffymetrixCelFile", function(this, unit=NULL) {
  .Deprecated("readUnits");
  this[units=unit, drop=TRUE];
}, protected=TRUE, deprecated=TRUE)


setMethodS3("[", "AffymetrixCelSet", function(this, units=NULL, ..., drop=FALSE) {
  .Deprecated("readUnits");
  res <- readUnits(this, units=units, ...);
  if (drop && length(res) == 1)
    res <- res[[1]];
  res;
}, protected=TRUE, deprecated=TRUE)

setMethodS3("[[", "AffymetrixCelSet", function(this, units=NULL, ...) {
  .Deprecated("readUnits");
  this[units=units, ..., drop=TRUE];
}, protected=TRUE, deprecated=TRUE)


setMethodS3("[", "CnagCfhFile", function(this, units=NULL, alleles=NULL, drop=FALSE) {
  .Deprecated("readUnits");
  data <- readUnits(this, units=units);
  if (!is.null(alleles)) {
    data <- data[, alleles, drop=drop];
  } else {
    if (drop && length(data) == 1)
      data <- data[[1]];
  }
  data;
}, protected=TRUE, deprecated=TRUE)

setMethodS3("[[", "CnagCfhFile", function(this, unit=NULL) {
  .Deprecated("readUnits");
  this[units=unit, drop=TRUE];
}, protected=TRUE, deprecated=TRUE)


setMethodS3("[", "CnagCfhSet", function(this, units=NULL, ..., drop=FALSE) {
  .Deprecated("readUnits");
  res <- readUnits(this, units=units, ...);
  if (drop && length(res) == 1)
    res <- res[[1]];
  res;
}, protected=TRUE, deprecated=TRUE)

setMethodS3("[[", "CnagCfhSet", function(this, units=NULL, ...) {
  .Deprecated("readUnits");
  this[units=units, ..., drop=TRUE];
}, protected=TRUE, deprecated=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fromName() => byName()
#
# 2009-09-05
# o CLEAN UP: Now static methods fromChipType() and fromName() of
#   AffymetrixCelSet and other classes are defunct.  Instead, use static
#   methods byChipType() and byName() instead.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("fromName", "AffymetrixCelSet", function(static, ...) {
  .Defunct("byName");
  className <- class(static)[1];
  msg <- sprintf("%s$fromName() is defunct. Use %s$byName() instead.",
                                                className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromName", "CnagCfhSet", function(static, ...) {
  .Defunct("byName");
  className <- class(static)[1];
  msg <- sprintf("%s$fromName() is defunct. Use %s$byName() instead.",
                                                className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromName", "DChipDcpSet", function(static, ...) {
  .Defunct("byName");
  className <- class(static)[1];
  msg <- sprintf("%s$fromName() is defunct. Use %s$byName() instead.",
                                                className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fromChipType() => byChipType()
#
# 2009-09-05
# o CLEAN UP: Now static methods fromChipType() and fromName() of
#   AffymetrixCelSet and other classes are defunct.  Instead, use static
#   methods byChipType() and byName() instead.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("fromChipType", "AffymetrixCsvGenomeInformation", function(static, ...) {
  .Defunct("byChipType");
  className <- class(static)[1];
  msg <- sprintf("%s$fromChipType() is defunct. Use %s$byChipType() instead.",
                                                        className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromChipType", "AffymetrixProbeTabFile", function(static, ...) {
  .Defunct("byChipType");
  className <- class(static)[1];
  msg <- sprintf("%s$fromChipType() is defunct. Use %s$byChipType() instead.",
                                                        className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromChipType", "AffymetrixTabularFile", function(static, ...) {
  .Defunct("byChipType");
  className <- class(static)[1];
  msg <- sprintf("%s$fromChipType() is defunct. Use %s$byChipType() instead.",
                                                        className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromChipType", "AffymetrixTsvFile", function(static, ...) {
  .Defunct("byChipType");
  className <- class(static)[1];
  msg <- sprintf("%s$fromChipType() is defunct. Use %s$byChipType() instead.",
                                                        className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromChipType", "AromaChipTypeAnnotationFile", function(static, ...) {
  .Defunct("byChipType");
  className <- class(static)[1];
  msg <- sprintf("%s$fromChipType() is defunct. Use %s$byChipType() instead.",
                                                        className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromChipType", "DChipGenomeInformation", function(static, ...) {
  .Defunct("byChipType");
  className <- class(static)[1];
  msg <- sprintf("%s$fromChipType() is defunct. Use %s$byChipType() instead.",
                                                        className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromChipType", "DChipSnpInformation", function(static, ...) {
  .Defunct("byChipType");
  className <- class(static)[1];
  msg <- sprintf("%s$fromChipType() is defunct. Use %s$byChipType() instead.",
                                                        className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromChipType", "GenomeInformation", function(static, ...) {
  .Defunct("byChipType");
  className <- class(static)[1];
  msg <- sprintf("%s$fromChipType() is defunct. Use %s$byChipType() instead.",
                                                        className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromChipType", "SnpInformation", function(static, ...) {
  .Defunct("byChipType");
  className <- class(static)[1];
  msg <- sprintf("%s$fromChipType() is defunct. Use %s$byChipType() instead.",
                                                        className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromChipType", "UflSnpInformation", function(static, ...) {
  .Defunct("byChipType");
  className <- class(static)[1];
  msg <- sprintf("%s$fromChipType() is defunct. Use %s$byChipType() instead.",
                                                        className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)

setMethodS3("fromChipType", "UgpGenomeInformation", function(static, ...) {
  .Defunct("byChipType");
  className <- class(static)[1];
  msg <- sprintf("%s$fromChipType() is defunct. Use %s$byChipType() instead.",
                                                        className, className);
  throw(msg);
}, static=TRUE, protected=TRUE, deprecated=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# readData() => readDataFrame()
#
# Version: 0.9.1.4 [2008-05-09] (Never released)
# o CLEAN UP: Renamed all readData() methods that return a data.frame
#   to readDataFrame().
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("readData", "AffymetrixTsvFile", function(this, ...) {
  .Defunct("readDataFrame");
  readDataFrame(this, ...);
}, protected=TRUE, deprecated=TRUE)

setMethodS3("readData", "GenomeInformation", function(this, ...) {
  .Defunct("readDataFrame");
  readDataFrame(this, ...);
}, protected=TRUE, deprecated=TRUE);

setMethodS3("readData", "SnpInformation", function(this, ...) {
  .Defunct("readDataFrame");
  readDataFrame(this, ...);
}, protected=TRUE, deprecated=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# getCellMap() => getUnitGroupCellMap()
#
# Version: 0.9.1.4 [2008-05-09] (Never released)
# o DEFUNCT: All getCellMap() are now defunct; use getUnitGroupCellMap().
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getCellMap", "ChipEffectFile", function(this, ...) {
  .Defunct("getUnitGroupCellMap");
  throw("getCellMap() is defunct. Use getUnitGroupCellMap() instead.");
}, protected=TRUE, deprecated=TRUE)

setMethodS3("getCellMap", "FirmaFile", function(this, ...) {
  .Defunct("getUnitGroupCellMap");
  throw("getCellMap() is defunct. Use getUnitGroupCellMap() instead.");
}, protected=TRUE, deprecated=TRUE)

setMethodS3("getCellMap", "ResidualFile", function(this, ...) {
  .Defunct("getUnitGroupCellMap");
  throw("getCellMap() is defunct. Use getUnitGroupCellMap() instead.");
}, protected=TRUE, deprecated=TRUE)

setMethodS3("getCellMap", "WeightsFile", function(this, ...) {
  .Defunct("getUnitGroupCellMap");
  throw("getCellMap() is defunct. Use getUnitGroupCellMap() instead.");
}, protected=TRUE, deprecated=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# getChipEffects() => getChipEffectSet()
# getProbeAffinities() => getProbeAffinityFile()
#
# Deprecated ca 2007.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getChipEffects", "ProbeLevelModel", function(this, ...) {
  .Defunct("getChipEffectSet");
  getChipEffectSet(this, ...);
}, protected=TRUE, deprecated=TRUE)

setMethodS3("getChipEffects", "QualityAssessmentModel", function(this, ...) {
  .Defunct("getChipEffectSet");
  getChipEffectSet(this, ...);
}, protected=TRUE, deprecated=TRUE)

setMethodS3("getChipEffects", "ExonRmaPlm", function(this, ...) {
  .Defunct("getChipEffectSet");
  getChipEffectSet(this, ...);
})

setMethodS3("getProbeAffinities", "ProbeLevelModel", function(this, ...) {
  .Defunct("getProbeAffinityFile");
  getProbeAffinityFile(this, ...);
}, protected=TRUE, deprecated=TRUE)

setMethodS3("getProbeAffinities", "ExonRmaPlm", function(this, ...) {
  .Defunct("getProbeAffinityFile");
  getProbeAffinityFile(this, ...);
}, protected=TRUE, deprecated=TRUE)

# 2008-09-03
# o Added getFitUnitGroupFunction() model, which is a better name than
#   getFitFunction().
setMethodS3("getFitFunction", "MultiArrayUnitModel", function(...) {
  .Defunct("getFitUnitGroupFunction");
  throw("getFitFunction() is deprecated. Please use getFitUnitGroupFunction() instead.");
}, protected=TRUE, deprecated=TRUE)


setMethodS3("getParameterSet", "Model", function(this, ...) {
  .Deprecated("getParameters");
  getParameters(this, ...);
}, protected=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Misc.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("bgAdjustOptical", "AffymetrixCelSet", function(this, ...) {
  .Defunct("OpticalBackgroundCorrection");
}, private=TRUE, deprecated=TRUE)

setMethodS3("bgAdjustRma", "AffymetrixCelSet", function(this, ...) {
  .Defunct("RmaBackgroundCorrection");
}, private=TRUE, deprecated=TRUE)


setMethodS3("getUnitSizes", "AffymetrixCdfFile", function(this, ...) {
  .Deprecated("nbrOfGroupsPerUnit");
  nbrOfGroupsPerUnit(this, ...);
}, protected=TRUE, deprecated=TRUE)


# 2013-04-29 [HB]
# o Made getUnique(), createUnique(), createMonoCell(), doCRMA() defunct.
# 2011-04-15 [HB]
# o DEPRECATED: getUnique() and createUnique() are deprecated.
#   Use getUniqueCdf() and createUniqueCdf() instead.
setMethodS3("getUnique", "AffymetrixCdfFile", function(this, ...) {
  .Defunct("getUniqueCdf");
  getUniqueCdf(this, ...);
}, protected=TRUE, deprecated=TRUE)


setMethodS3("createUnique", "AffymetrixCdfFile", function(this, ...) {
  .Defunct("createUniqueCdf");
  createUniqueCdf(this, ...);
}, protected=TRUE, deprecated=TRUE)


setMethodS3("getMonoCell", "AffymetrixCdfFile", function(this, ...) {
  .Deprecated("getMonocellCdf");
  getMonocellCdf(this, ...);
}, protected=TRUE, deprecated=TRUE)


setMethodS3("createMonoCell", "AffymetrixCdfFile", function(this, ...) {
  .Defunct("createMonocellCdf");
  createMonocellCdf(this, ...);
}, protected=TRUE, deprecated=TRUE)


setMethodS3("getMatrixChipEffectFiles", "CopyNumberChromosomalModel", function(...) {
  .Deprecated("getDataFileMatrix");
  getDataFileMatrix(...);
}, protected=TRUE, deprecated=TRUE)


setMethodS3("calculateResiduals", "ProbeLevelModel", function(this, ...) {
  .Deprecated("calculateResidualSet");
  calculateResidualSet(this, ...);
}, private=TRUE)


setMethodS3("calculateResiduals", "FirmaModel", function(this, ...) {
  .Deprecated("calculateResidualSet");
  calculateResidualSet(this, ...);
}, private=TRUE)


# 2011-11-19
# o Deprecated doCRMA().  Use doCRMAv1() or doCRMAv2() instead.
setMethodS3("doCRMA", "default", function(dataSet, chipTypes=NULL, ..., logName=NULL, ram=NULL, verbose=-8) {
  .Defunct("doCRMAv1");

  csRawList <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (inherits(dataSet, "AffymetrixCelSet")) {
    csRawList <- list(dataSet);
    chipTypes <- sapply(csRawList, FUN=function(cs) {
      getChipType(getCdf(cs));
    });
    names(csRawList) <- chipTypes;
    dataSet <- getFullName(csRawList[[1]]);
  } else if (is.list(dataSet)) {
    csRawList <- dataSet;
    # Validate elements
    lapply(csRawList, FUN=function(cs) {
      cs <- Arguments$getInstanceOf(cs, "AffymetrixCelSet", .name="dataSet");
    });
    chipTypes <- sapply(csRawList, FUN=function(cs) {
      getChipType(getCdf(cs));
    });
    names(csRawList) <- chipTypes;
    dataSet <- getFullName(csRawList[[1]]);
  } else {
    dataSet <- Arguments$getCharacter(dataSet, length=c(1,1));
  }

  # Argument 'chipTypes':
  if (!is.null(chipTypes)) {
    chipTypes <- Arguments$getCharacters(chipTypes);
  }

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Return object
  fit <- list();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup log file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(logName)) {
    mergedChipTypes <- mergeByCommonTails(chipTypes, collapse="+");
    logName <- sprintf("%s,%s", dataSet, mergedChipTypes);
  }
  logDate <- format(Sys.time(), "%Y%m%d-%H%M%S");
  logPathname <- sprintf("%s,%s.log", logName, logDate);
  verbose && cat(verbose, "Log file: ", logPathname);

  fileLog <- Verbose(con=logPathname, threshold=-50, timestamp=TRUE);
  fileLog && header(fileLog, "Log file created by doCRMA() in the aroma.affymetrix package");


  # Setup multi-verbose output
  log <- MultiVerbose(list(fileLog, verbose), threshold=-100);
  timestampOn(log);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  log && header(log, "Setting up raw data set(s)");
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(csRawList)) {
    csRawList <- list();
    for (chipType in chipTypes) {
      cs <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
      log && print(log, cs);
      csRawList[[chipType]] <- cs;
      rm(cs);
    }
  }
  log && print(log, csRawList);
  fit$csRawList <- csRawList;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  log && header(log, "Asserting that all annotation data needed is available");
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (chipType in chipTypes) {
    log && enter(log, "Chip type: ", chipType);

    # Is chip type supported?
    if (regexpr("^Mapping", chipType) != -1) {
    } else if (regexpr("^GenomeWideSNP_", chipType) != -1) {
    } else {
      throw("Unsupported chip type");
    }

    cs <- csRawList[[chipType]];
    log && print(log, cs);

    # Get the CDF
    cdf <- getCdf(cs);
    log && print(log, cdf);

    # Get the genome information file
    gi <- getGenomeInformation(cdf);
    log && print(log, gi);

    # Get the SNP information file
    si <- getSnpInformation(cdf);
    log && print(log, si);

    log && exit(log);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  log && header(log, "Calibrating for allelic crosstalk");
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  csCList <- list();
  for (chipType in chipTypes) {
    log && enter(log, "Chip type: ", chipType);

    # Input data
    cs <- csRawList[[chipType]];

    # Annotation data
    cdf <- getCdf(cs);
    gi <- getGenomeInformation(cdf);

    # Identify units on sex chromosomes
    chromosomes <- intersect(23:24, getChromosomes(gi));
    unitsXY <- getUnitsOnChromosomes(gi, chromosomes);
    nbrOfUnits <- nbrOfUnits(cdf);
    log && printf(log, "Identified %d units (%.1f%%) on sex chromosomes out of %d\n", length(unitsXY), 100*length(unitsXY)/nbrOfUnits, nbrOfUnits);

    cellsXY <- getCellIndices(cdf, units=unitsXY, useNames=FALSE, unlist=TRUE);
    rm(unitsXY);
    nbrOfCells <- nbrOfCells(cdf);
    cellsNotXY <- setdiff(1:nbrOfCells, cellsXY);
    log && printf(log, "Identified %d cells (%.1f%%) on sex chromosomes out of %d\n", length(cellsXY), 100*length(cellsXY)/nbrOfCells, nbrOfCells);
    rm(cellsXY);

    tags <- c("*", "-XY");
    acc <- AllelicCrosstalkCalibration(cs, subsetToAvg=cellsNotXY, tags=tags);
    rm(cellsNotXY);
    log && print(log, acc);

    # Output data
    csC <- process(acc, verbose=log);
    log && print(log, csC);

    # Validation
    stopifnot(identical(getNames(csC), getNames(cs)));

    # Storing results
    csCList[[chipType]] <- csC;
    rm(cs, csC);

    # Garbage collect
    gc <- gc();
    log && print(log, gc);

    log && exit(log);
  }
  fit$csCList <- csCList;



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  log && header(log, "Summarizing probe signals");
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cesList <- list();
  for (chipType in chipTypes) {
    log && enter(log, "Chip type: ", chipType);

    # Input data
    csC <- csCList[[chipType]];

    # Setting up the probe-level model
    plmModel <- NULL;
    if (regexpr("^Mapping", chipType) != -1) {
      plmModel <- RmaCnPlm;
    } else if (regexpr("^GenomeWideSNP_", chipType) != -1) {
      plmModel <- AvgCnPlm;
    }

    plmShift <- +300;
    plmShift <- as.integer(plmShift);
    tags <- c(sprintf("%+d", plmShift), "*");
    plm <- plmModel(csC, mergeStrands=TRUE, combineAlleles=TRUE, tags=tags);
    plm$shift <- as.double(plmShift);
    plm$treatNAsAs <- "weights";
    log && print(log, plm);

    # Fitting probe-level model
    fit(plm, ram=ram, verbose=log);

    # Getting the chip effects
    ces <- getChipEffectSet(plm);
    log && print(log, ces);
    rm(plm);

    # Validation
    stopifnot(identical(getNames(ces), getNames(csC)));

    # Storing results
    cesList[[chipType]] <- ces;
    rm(csC, ces);

    # Garbage collect
    gc <- gc();
    log && print(log, gc);

    log && exit(log);
  }
  fit$cesList <- cesList;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  log && header(log, "Normalizing for fragment-length effects");
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cesNList <- list();
  for (chipType in chipTypes) {
    log && enter(log, "Chip type: ", chipType);

    # Input data
    ces <- cesList[[chipType]];

    # Annotation data
    cdf <- getCdf(ces);
    gi <- getGenomeInformation(cdf);

    # Identify units on sex chromosomes
    chromosomes <- intersect(23:24, getChromosomes(gi));
    unitsXY <- getUnitsOnChromosomes(gi, chromosomes);
    nbrOfUnits <- nbrOfUnits(cdf);
    unitsNotXY <- setdiff(1:nbrOfUnits, unitsXY);
    log && printf(log, "Identified %d units (%.1f%%) on sex chromosomes out of %d\n", length(unitsXY), 100*length(unitsXY)/nbrOfUnits, nbrOfUnits);
    rm(unitsXY);

    # Setting up the fragment-length normalization model
    tags <- c("*", "-XY");
    fln <- FragmentLengthNormalization(ces, subsetToFit=unitsNotXY, tags=tags);
    log && print(log, fln);
    rm(unitsNotXY);

    # Processing
    cesN <- process(fln, verbose=log);
    log && print(log, cesN);
    rm(fln);

    # Validation
    stopifnot(identical(getNames(cesN), getNames(ces)));

    # Storing results
    cesNList[[chipType]] <- cesN;
    rm(ces, cesN);

    # Garbage collect
    gc <- gc();
    log && print(log, gc);

    log && exit(log);
  }
  fit$cesNList <- cesNList;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  log && header(log, "Returning results");
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit;
}, protected=TRUE, deprecated=TRUE)



setMethodS3("getExpectedOutputFiles", "MatSmoothing", function(this, ...) {
  .Deprecated("getExpectedOutputFullnames");
  fullnames <- getExpectedOutputFullnames(this, ...);

  # "Dummy" filenames
  filenames <- sprintf("%s.CEL", fullnames);

  filenames;
}, protected=TRUE, deprecated=TRUE)


# 2013-04-29
# o Made PdInfo2Cdf() defunct.
# 2012-03-23
# o CLEANUP: Deprecated PdInfo2Cdf() in favor (identical) pdInfo2Cdf(),
#   because the former does not follow the Aroma naming conventions.
PdInfo2Cdf <- function(...) {
  .Defunct("pdInfo2Cdf");
  throw("PdInfo2Cdf() is deprecated. Use pdInfo2Cdf() instead, which works identically.");
  pdInfo2Cdf(...);
} # PdInfo2Cdf()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TODO, but still used alot internally.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getData", "AffymetrixCelFile", function(this, ...) {
##  .Deprecated("readRawData");
  readRawData(this, ...);
}, protected=TRUE, deprecated=TRUE)



############################################################################
# HISTORY:
# 2013-04-29
# o CLEANUP: Made several deprecated methods defunct.
# 2012-11-20
# o Deprecated getParameterSet() in favor of old getParameters().
# o Defuncted bgAdjust(Optical|Rma)() for AffymetrixCelSet.
# 2012-10-17
# o Deprecated getMonoCell() and createMonoCell().
# 2012-10-14
# o Created 999.DEPRECATED.R.
# 2011-02-19
# o Replaced deprecated getListOfChipEffectSets() with getSets() for
#   ChromosomalModel and SmoothMultiarrayModel.
# 2009-09-05
# o CLEAN UP: Now static methods fromChipType() and fromName() of
#   AffymetrixCelSet and other classes are defunct.  Instead, use static
#   methods byChipType() and byName() instead.
# 2008-06-06
# o Removed deprecated getTargetFunction() for FragmentLengthNormalization.
############################################################################
