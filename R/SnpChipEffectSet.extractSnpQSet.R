setMethodS3("extractSnpQSet", "SnpChipEffectSet", function(this, units=NULL, sortUnits=TRUE, transform=log2, ..., verbose=FALSE) {
  require("oligo") || throw("Package not loaded: oligo");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this, fullname=FALSE);

  if (regexpr("^GenomeWideSNP_(5|6)$", chipType) != -1) {
    throw("Cannot extract SnpQSet: Unsupported chip type: ", chipType);
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
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
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
#  dirs <- lapply(dirs, FUN=function(groups) groups[1]);
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
  # Allocate and populate SnpQSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocate and populate SnpQSet");
  res <- new("SnpQSet", 
    senseThetaA     = transform(theta[,1,,drop=TRUE]),
    senseThetaB     = transform(theta[,2,,drop=TRUE]),
    antisenseThetaA = transform(theta[,3,,drop=TRUE]),
    antisenseThetaB = transform(theta[,4,,drop=TRUE])
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
# 2008-12-05
# o Created.
############################################################################
