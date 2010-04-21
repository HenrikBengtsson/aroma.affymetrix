###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod groupUnitsByDimension
#
# @title "Groups units by dimensions"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{An optional @integer @vector specifying the units to be
#     queried.}
#   \item{...}{Not used.}
#   \item{sort}{If @TRUE, the sets are ordered by number of groups per
#     units and then by the number of cells per group.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a nested @list.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("groupUnitsByDimension", "AffymetrixCdfFile", function(this, units=NULL, ..., sort=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
    allUnits <- units;
  } else {
    allUnits <- seq(length=nbrOfUnits(this));
  }

  # Argument 'sort':
  sort <- Arguments$getLogical(sort);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Grouping units by dimensions");

  verbose && cat(verbose, "Units: ");
  verbose && str(verbose, allUnits);

  # Identify SNPs with equal number of cells
  sizes <- nbrOfCellsPerUnitGroup(this, units=units, verbose=less(verbose,1));

  # Identify sets of units of equal dimensions
  nbrOfGroupsPerUnit <- sapply(sizes, FUN=length);
  t <- table(nbrOfGroupsPerUnit);
  verbose && print(verbose, t);
  ## nbrOfGroupsPerUnit
  ##      2      4
  ## 906874   2748

  uSizes <- as.integer(names(t));

  # Sort?
  if (sort) {
    uSizes <- sort(uSizes);
  }

  res <- list();

  for (uu in seq(along=uSizes)) {
    sizeUU <- uSizes[uu];
    verbose && enter(verbose, sprintf("Unit size %d of %d", uu, length(uSizes)));
    verbose && cat(verbose, "Number of groups per unit: ", sizeUU);

    # Pull out all units 'size' number of groups
    idxsUU <- which(nbrOfGroupsPerUnit == sizeUU);
    unitsUU <- allUnits[idxsUU];
    verbose && cat(verbose, "Units: ");
    verbose && str(verbose, unitsUU);

    # Group all those units into those with equal number cells per group
    sizesUU <- sizes[idxsUU];
    sizesUU <- unlist(sizesUU, use.names=FALSE);
    sizesUU <- matrix(sizesUU, ncol=sizeUU, byrow=TRUE);
    uDimsUU <- unique(sizesUU);

    # Sort?
    if (sort) {
      o <- order(uDimsUU[,1]);
      uDimsUU <- uDimsUU[o,,drop=FALSE];
    }

    verbose && cat(verbose, "Unique dimensions: ");
    verbose && print(verbose, uDimsUU);

    resUU <- list();
    resUU$nbrOfGroups <- sizeUU;
    resUU$units <- unitsUU;
    resUU$sets <- list();

    # For each group of equal dimension
    for (kk in seq(length=nrow(uDimsUU))) {
      dimKK <- uDimsUU[kk, ,drop=TRUE];
      verbose && enter(verbose, sprintf("Dimension %d of %d", kk, nrow(uDimsUU)));
      verbose && cat(verbose, "Number of cells per group: ", paste(dimKK, collapse=", "));

      keep <- rep(TRUE, times=length(unitsUU));
      for (cc in seq(length=ncol(uDimsUU))) {
        keep <- keep & (dimKK[cc] == sizesUU[,cc, drop=TRUE]);
      }
      idxsKK <- which(keep);
      unitsKK <- unitsUU[idxsKK];
#      verbose && cat(verbose, "Units: ");
#      verbose && str(verbose, unitsKK);

      resKK <- list();
      resKK$nbrOfCellsPerGroup <- dimKK;
      resKK$units <- unitsKK;
      verbose && str(verbose, resKK);

      resUU$sets[[kk]] <- resKK;
      rm(idxsKK, unitsKK, resKK);

      verbose && exit(verbose);
    } # for (kk ...)

    res[[uu]] <- resUU;
    rm(keep, idxsUU, unitsUU, resUU);

    verbose && exit(verbose);
  } # for (uu ...) 

  verbose && exit(verbose);

  res;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2010-04-21
# o Created.
############################################################################
