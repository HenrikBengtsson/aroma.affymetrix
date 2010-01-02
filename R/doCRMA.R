setMethodS3("doCRMA", "default", function(dataSet, chipTypes=NULL, ..., logName=NULL, ram=NULL, verbose=-8) {
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
})


############################################################################
# HISTORY:
# 2007-11-25
# o Created.
############################################################################
