setMethodS3("extractTheta", "ChipEffectSet", function(this, units=NULL, groups=NULL, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdf <- getCdf(this);

  # Argument 'units':
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    ugcMap <- NULL;
  } else if (isUnitGroupCellMap(units)) {
    ugcMap <- units;
    units <- unique(ugcMap[,"unit"]);
    nbrOfUnits <- length(units);
    rm(units);
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
    nbrOfUnits <- length(units);
    ugcMap <- NULL;
  }

  # Argument 'groups':
  if (!is.null(groups)) {
    groups <- Arguments$getIndices(groups, max=999);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the UGC map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }

  if (!is.null(groups)) {
    idxs <- which(ugcMap$group %in% groups);
    ugcMap <- ugcMap[idxs,,drop=FALSE];
  } else {
    groups <- sort(unique(ugcMap$group));
  }

  verbose && cat(verbose, "Using (unit,group,cell) map:");
  verbose && str(verbose, ugcMap);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the thetas
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  nbrOfGroups <- length(groups);
  nbrOfArrays <- nbrOfArrays(this);
  dim <- c(nbrOfUnits, nbrOfGroups, nbrOfArrays);
  dimnames <- list(NULL, NULL, getNames(this));
  naValue <- as.double(NA);
  theta <- array(naValue, dim=dim, dimnames=dimnames);
  for (kk in seq(length=nbrOfArrays)) {
    ce <- getFile(this, kk);
    thetaKK <- extractTheta(ce, units=ugcMap, groups=groups, 
                                                 verbose=less(verbose, 5));
    verbose && str(verbose, thetaKK);
    theta[,,kk] <- thetaKK;
  }
  rm(ugcMap);

  # Drop singleton dimensions
  if (drop) {
    theta <- drop(theta);
  }

  verbose && cat(verbose, "Thetas:");
  verbose && str(verbose, theta);

  theta;
})



setMethodS3("extractTheta", "SnpChipEffectSet", function(this, groups=NULL, ...) {
  if (is.null(groups)) {
    maxNbrOfGroups <- 4;
    if (this$mergeStrands) {
      maxNbrOfGroups <- maxNbrOfGroups / 2;
    }
    groups <- 1:maxNbrOfGroups;
  }

  theta <- NextMethod("extractTheta", this, groups=groups, ...);

  theta;
})



setMethodS3("extractTheta", "CnChipEffectSet", function(this, groups=NULL, ...) {
  if (is.null(groups)) {
    maxNbrOfGroups <- 4;
    if (this$mergeStrands) {
      maxNbrOfGroups <- maxNbrOfGroups / 2;
    }
    if (this$combineAlleles) {
      maxNbrOfGroups <- maxNbrOfGroups / 2;
    }
    groups <- 1:maxNbrOfGroups;
  }

  theta <- NextMethod("extractTheta", this, groups=groups, ...);

  theta;
})


############################################################################
# HISTORY:
# 2008-07-13
# o Added argument 'drop=FALSE' to extractTheta().
# 2008-05-10
# o Created.
############################################################################
