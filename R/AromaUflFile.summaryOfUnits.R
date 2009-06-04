setMethodS3("summaryOfUnits", "AromaUflFile", function(this, enzymeLabels=paste("enzyme", 1:nbrOfEnzymes(this), sep=""), unitClasses=c(snp="^SNP_", cnp="^CN_", affxSnp="^AFFX-SNP_"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  unf <- getUnitNamesFile(this);
  nbrOfEnzymes <- nbrOfEnzymes(this);

  # Argument 'enzymeLabels':
  enzymeLabels <- Arguments$getCharacters(enzymeLabels, length=nbrOfEnzymes);

  # Argument 'unitClasses':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Summarizing UFL data by unit and enzyme classes");

  verbose && enter(verbose, "Extracting fragment-length information");
  # Extract fragment-length data
  fl <- as.matrix(this[]);
  colnames(fl) <- enzymeLabels;
  verbose && exit(verbose);

  # Check for existing data
  hasFl <- is.finite(fl);

  verbose && enter(verbose, "Identifying unit classes");
  # Extract unit classes of interest
  patterns <- unitClasses;
  unitClasses <- lapply(patterns, FUN=function(pattern) {
    indexOf(unf, pattern);
  })
  names(unitClasses) <- names(patterns);
  nbrOfUnits <- nbrOfUnits(unf);
  unitClasses[["other"]] <- setdiff(1:nbrOfUnits, 
                              unlist(unitClasses, use.names=FALSE));
  verbose && exit(verbose);
 
  verbose && enter(verbose, "Identifying enzyme classes");
  # Extract enzyme classes of interest
  enzymeClasses <- list();
  ees <- 1:ncol(hasFl);
  # Identify units exclusively on one enzyme
  for (ee in ees) {
    ok <- hasFl[,ee];
    for (ff in setdiff(ees, ee)) {
      ok <- ok & (!hasFl[,ff]);
    }
    enzymeClasses[[ee]] <- which(ok);
  }
  names(enzymeClasses) <- paste(enzymeLabels, "-only", sep="");

  # Identify units that are all enzymes
  if (nbrOfEnzymes > 1) {
    ok <- rep(TRUE, nrow(hasFl));
    for (ee in ees)
      ok <- ok & (hasFl[,ee]);
    name <- ifelse(nbrOfEnzymes == 2, "both", "all")
    enzymeClasses[[name]] <- which(ok);
  }

  # Identifying units for which there is no data
  nas <- rep(TRUE, nrow(hasFl));
  for (ee in ees)
    nas <- nas & (!hasFl[,ee]);
  enzymeClasses[["missing"]] <- which(nas);
  verbose && exit(verbose);


  # Build (unit,enzyme) sets of interest
  verbose && enter(verbose, "Building (unit, enzyme) classes");
  unitSets <- vector("list", length(unitClasses))
  names(unitSets) <- names(unitClasses);
  for (ii in seq(along=unitClasses)) {
    unitSet <- vector("list", length(enzymeClasses))
    names(unitSet) <- names(enzymeClasses)
    for (jj in seq(along=enzymeClasses)) {
      units <- intersect(unitClasses[[ii]], enzymeClasses[[jj]])
      unitSet[[jj]] <- units
    }
    unitSets[[ii]] <- unitSet
  }
  verbose && exit(verbose);
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Summary table of different unit types
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Tabulating data");
  utbl <- sapply(unitSets, FUN=sapply, length)
  utbl <- cbind(utbl, total=rowSums(utbl))
  utbl <- rbind(utbl, total=colSums(utbl))
  verbose && exit(verbose);

  verbose && exit(verbose);

  utbl;
}) # summaryOfUnits()



############################################################################
# HISTORY:
# 2009-05-20
# o Updated to make use of getUnitNamesFile() instead of getCdf().
# 2007-12-17
# o Created.
############################################################################

