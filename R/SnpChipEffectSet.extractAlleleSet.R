setMethodS3("extractAlleleSet", "SnpChipEffectSet", function(this, units=NULL, sortUnits=TRUE, transform=log2, ..., verbose=FALSE) {
  require("Biobase") || throw("Package not loaded: Biobase");
  require("oligo") || throw("Package not loaded: oligo");

  # Assert oligo version
  pkg <- Package("oligo");
  if (isOlderThan(Package("oligo"), "1.12.0")) {
    throw("extractAlleleSet() requires oligo v1.12.0 or newer: ", getVersion(pkg));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this, fullname=FALSE);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'units':  
  if (is.null(units)) {
    # Identify all SNP_A-* units (as what is returned by oligo)
    units <- indexOf(cdf, pattern="^SNP_A-");
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Order units lexicographically by their names?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  unitNames <- getUnitNames(cdf, units=units);
  if (sortUnits) {
    verbose && enter(verbose, "Sort units by their names");
    srt <- sort(unitNames, method="quick", index.return=TRUE);
    unitNames <- srt$x;
    units <- units[srt$ix];
    rm(srt);  # Not needed anymore
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Inferring what unit groups to extract
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Inferring the number of groups to extract");
  verbose && cat(verbose, "Chip type: ", getChipType(this));
  ugcMap <- getUnitGroupCellMap(this, units=units, verbose=verbose);
  maxGroup <- max(ugcMap$group, na.rm=TRUE);
  rm(ugcMap);
  verbose && cat(verbose, "Max number of groups unit (unit,group,cell) map: ", maxGroup);
  if (maxGroup > 4) {
    maxGroup <- 4L;
  }

  # Sanity check
  if (!is.element(maxGroup, c(2,4))) {
    throw("Unsupported value on 'maxGroup': ", maxGroup);
  }

  groups <- 1:maxGroup;

  verbose && cat(verbose, "Groups to extract: ", seqToHumanReadable(groups));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Infer which (antisense, sense) unit groups to swap (more below)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (maxGroup == 4) {
    verbose && enter(verbose, "Inferring which unit groups to be swapped to (sense, antisense)");
  
    # Make sure pairs are order as (sense, antisense)
    dirs <- getGroupDirections(cdf, units=units, verbose=less(verbose, 5));
    names(dirs) <- NULL;
  
    gc <- gc();
    verbose && print(verbose, gc);
  
    # Sanity check
    lens <- sapply(dirs, FUN=length);
    uLens <- unique(lens);
    if (any(!is.element(uLens, unique(c(2,maxGroup))))) {
      throw("Internal error: Unexpected number of unit groups: ", 
                                                paste(uLens, collapse=", "));
    }
  
    # Extract the direction/strand of the first group
    # dirs <- lapply(dirs, FUN=function(groups) groups[1]);
    dirs <- lapply(dirs, FUN=.subset, 1);
    dirs <- unlist(dirs, use.names=FALSE);
  
    # Identify which to swap from (antisense,sense) to (sense,antisense)
    idxs <- which(dirs == 2);
    rm(dirs);

    gc <- gc();
    verbose && print(verbose, gc);
  
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting data");
  theta <- extractTheta(this, groups=groups, units=units, verbose=verbose);

  # Not needed anymore
  rm(units);

  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Swap (antisense, sense) unit groups to (sense, antisense)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (maxGroup == 4 && length(idxs) > 0) {
    verbose && enter(verbose, "Ordering unit groups to be (sense, antisense)");
    verbose && cat(verbose, "Swapping elements:");
    verbose && str(verbose, idxs);
    theta[idxs,,] <- theta[idxs,c(3,4,1,2),,drop=FALSE];
    rm(idxs);   # Not needed anymore
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Transform data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(transform)) {
    verbose && enter(verbose, "Transforming signals");
    theta <- transform(theta);
    verbose && str(verbose, theta);
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate and populate AlleleSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocate and populate AlleleSet");
  if (maxGroup == 2) {
    res <- new("AlleleSet", 
      alleleA = theta[,1,,drop=TRUE],
      alleleB = theta[,2,,drop=TRUE]
    );
  } else if (maxGroup == 4) {
    res <- new("AlleleSet", 
      antisenseAlleleA = theta[,3,,drop=TRUE],
      senseAlleleA     = theta[,1,,drop=TRUE],
      antisenseAlleleB = theta[,4,,drop=TRUE],
      senseAlleleB     = theta[,2,,drop=TRUE]
    );
  }

  # Not needed anymore
  rm(theta);

  # Assign feature data
  featureNames(res) <- unitNames;
  rm(unitNames);

  # Assign annotation data
  pdPkgName <- oligo::cleanPlatformName(chipType);
  annotation(res) <- pdPkgName;

  # Assign sample names
  filenames <- sapply(this, getFilename);
  names(filenames) <- NULL;
  filenames <- gsub(",chipEffects", "", filenames);
  sampleNames(res) <- filenames;

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2012-09-01
# o ROBUSTNESS: extractAlleleSet() for SnpChipEffectSet would throw an
#   exception if the 'Biobase' package was not loaded.
# 2010-05-11
# o BUG FIX: Too early rm(units) in extractAlleleSet().
# 2010-05-09
# o Now extractAlleleSet() handles SNP platforms with 2 & 4 unit groups.
# 2010-05-06
# o extractAlleleSet() now asserts that oligo v1.12.0 or older is installed.
# o Created from extractSnpQSet().
############################################################################
