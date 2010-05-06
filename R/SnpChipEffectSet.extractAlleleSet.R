setMethodS3("extractAlleleSet", "SnpChipEffectSet", function(this, units=NULL, sortUnits=TRUE, transform=log2, ..., verbose=FALSE) {
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

  # TO DO: Relax this test, because it should work. /HB 2010-05-06
  if (regexpr("^GenomeWideSNP_(5|6)$", chipType) != -1) {
    throw("Cannot extract AlleleSet: Unsupported chip type: ", chipType);
  }

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
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting data");
  theta <- extractTheta(this, groups=1:4, units=units, verbose=verbose);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Swap (antisense, sense) unit groups to (sense, antisense)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Ordering unit groups to be (sense, antisense)");
  # Make sure pairs are order as (sense, antisense)
  dirs <- getGroupDirections(cdf, units=units);
  names(dirs) <- NULL;
  rm(units);

  # Sanity check
  lens <- sapply(dirs, FUN=length);
  uLens <- unique(lens);
  if (any(!is.element(uLens, c(2,4)))) {
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

  verbose && cat(verbose, "Swapping elements:");
  verbose && str(verbose, idxs);
  theta[idxs,,] <- theta[idxs,c(3,4,1,2),];
  rm(idxs);   # Not needed anymore
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Transform data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(transform)) {
    verbose && enter(verbose, "Transforming signals");
    theta <- transform(theta);
    verbose && str(verbose, theta);
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate and populate AlleleSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocate and populate AlleleSet");
  res <- new("AlleleSet", 
    antisenseAlleleA = theta[,3,,drop=TRUE],
    senseAlleleA     = theta[,1,,drop=TRUE],
    antisenseAlleleB = theta[,4,,drop=TRUE],
    senseAlleleB     = theta[,2,,drop=TRUE]
  );

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
# 2010-05-06
# o extractAlleleSet() now asserts that oligo v1.12.0 or older is installed.
# o Created from extractSnpQSet().
############################################################################
