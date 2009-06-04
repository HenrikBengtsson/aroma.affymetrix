###########################################################################/**
# @RdocClass ExonChipEffectFile
#
# @title "The ExonChipEffectFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in the probe-level models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ChipEffectFile".}
#   \item{mergeGroups}{Specifies if the groups are merged or not for these
#      estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#   An object of this class is typically part of a @see "ExonChipEffectSet".
# }
#*/###########################################################################
setConstructorS3("ExonChipEffectFile", function(..., mergeGroups=FALSE) {
  this <- extend(ChipEffectFile(...), "ExonChipEffectFile",
    "cached:.cellIndices" = NULL,
    mergeGroups = mergeGroups
  );

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})


setMethodS3("clearCache", "ExonChipEffectFile", function(this, ...) {
  # Clear all cached values.
  for (ff in c(".cellIndices")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)



setMethodS3("getParameters", "ExonChipEffectFile", function(this, ...) {
  params <- NextMethod(generic="getParameters", object=this, ...);
  params$mergeGroups <- this$mergeGroups;
  params;
})


setMethodS3("getCellIndices", "ExonChipEffectFile", function(this, ..., force=FALSE, .cache=TRUE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "getCellIndices.ExonChipEffectFile()");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force || .cache) {
    chipType <- getChipType(getCdf(this));
    params <- getParameters(this);
    key <- list(method="getCellIndices", class=class(this)[1], 
                pathname=getPathname(this),
                chipType=chipType, params=params, ...);
    dirs <- c("aroma.affymetrix", chipType);
    id <- digest2(key);
  }

  if (!force) {
    # In memory?
    res <- this$.cellIndices[[id]];
    # On file?
    if (is.null(res)) {
      res <- loadCache(key=list(id), dirs=dirs);
      if (!is.null(res))
        where <- "on file";
    } else {
      where <- "in memory";
    }
    if (!is.null(res)) {
      size <- object.size(res);
      verbose && printf(verbose, "Returning value cached %s: %.1fMB\n", 
                                                   where, size/1024^2);
      verbose && exit(verbose);
      return(res);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get and restructure cell indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cells <- NextMethod("getCellIndices", this, ..., force=force, 
                                          .cache=FALSE, verbose=verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Merge groups?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If merging groups, we only need one chip-effect parameter per unit
  if (this$mergeGroups) {
    verbose && enter(verbose, "Merging groups");
    cells <- applyCdfGroups(cells, function(groups) {
      .subset(groups, 1);
    })
    verbose && exit(verbose);
  }

  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.cache) {
    # In-memory or on-file cache?
    if (object.size(cells) < 10e6) { 
      # In-memory cache for objects < 10Mb.
      this$.cellIndices <- list();
      this$.cellIndices[[id]] <- cells;
      verbose && cat(verbose, "Cached in memory");
    } else {
      # On-file cache
      # Keep, in-memory cache.
      if (!is.list(this$.cellIndices))
        this$.cellIndices <- list();
      this$.cellIndices[[id]] <- NULL;
      saveCache(cells, key=list(id), dirs=dirs);
      verbose && cat(verbose, "Cached to file");
    }
  }

  verbose && exit(verbose);

  cells;
})


setMethodS3("readUnits", "ExonChipEffectFile", function(this, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Check for cached data
  key <- list(method="readUnits", class=class(this)[1], 
              mergeGroups=this$mergeGroups, ...);
  id <- digest2(key);
  res <- this$.readUnitsCache[[id]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "readUnits.ExonChipEffectFile(): Returning cached data");
    return(res);
  }

  # Retrieve the data
  res <- NextMethod("readUnits", this, ..., force=TRUE, cache=FALSE, verbose=verbose);

  # Store read units in cache?
  if (cache) {
    verbose && cat(verbose, "readUnits.ExonChipEffectFile(): Updating cache");
    this$.readUnitsCache <- list();
    this$.readUnitsCache[[id]] <- res;
  }

  res;
})


## setMethodS3("findUnitsTodo", "ExonChipEffectFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
##   # Argument 'verbose':
##   verbose <- Arguments$getVerbose(verbose);
##   if (verbose) {
##     pushState(verbose);
##     on.exit(popState(verbose));
##   }
## 
## 
##   verbose && enter(verbose, "Identifying non-fitted units in chip-effect file");
##   verbose && cat(verbose, "Pathname: ", getPathname(this));
## 
## 
##   idxs <- NULL;
##   if (is.null(units)) {
##     # Look up chip-type and parameter specific but data set independent data
##     cdf <- getCdf(this);
##     chipType <- getChipType(cdf);
##     key <- list(method="findUnitsTodo", class=class(this)[1],
##                 chipType=chipType, params=getParameters(this));
##     dirs <- c("aroma.affymetrix", chipType);
##     if (!force) {
##       idxs <- loadCache(key, dirs=dirs);
##       if (!is.null(idxs))
##         verbose && cat(verbose, "Found indices cached on file");
##     }
##   }
## 
##   if (is.null(idxs)) {
##     verbose && enter(verbose, "Identifying CDF units");
## 
##     verbose && enter(verbose, "Reading CDF cell indices");
##     idxs <- getCellIndices(this, units=units, verbose=less(verbose));
##     # Example:
##     #  $ 2315554:List of 1
##     #   ..$ groups:List of 1
##     #   .. ..$ groupA:List of 1
##     #   .. .. ..$ indices: int 1
##     #   .. ..$ groupB:List of 1
##     #   .. .. ..$ indices: int 23
##     verbose && exit(verbose);
## 
##     verbose && enter(verbose, "Extracting first CDF block for each unit");
##     idxs <- applyCdfGroups(idxs, .subset2, 1);
##     # Example:
##     #  $ 2315554:List of 1
##     #   ..$ groups:List of 1
##     #   .. ..$ indices: int 1
##     ## The below makes no difference.  Thus, this function give exactly
##     ## the same result as the one in the super class. Note, mergeGroups
##     ## has already been taken care of inside getCellIndices() /HB 2007-08-17
##     if (this$mergeGroups) {
##       # do a further round of reduction
##       idxs <- applyCdfGroups(idxs, .subset2, 1);
##       # Example:
##       #  $ 2315554:List of 1
##       #   ..$ groups: int 1
##     }
##     verbose && exit(verbose);
## 
##     idxs <- unlist(idxs, use.names=FALSE);
## 
##     if (is.null(units)) {
##       verbose && enter(verbose, "Saving to file cache");
##       saveCache(idxs, key=key, dirs=dirs);
##       verbose && exit(verbose);
##     }
## 
##     verbose && exit(verbose);
##   }
## 
## 
##   # Read one cell from each unit
##   verbose && enter(verbose, "Reading data for these ", length(idxs), " cells");
##   value <- readCel(getPathname(this), indices=idxs, readIntensities=FALSE,
##                    readStdvs=TRUE, readPixels=FALSE)$stdvs;
##   verbose && exit(verbose);
## 
## 
##   # Identify units for which the stdvs <= 0.
##   value <- which(value <= 0);
## 
##   if (!is.null(units))
##     value <- units[value];
##   verbose && cat(verbose, "Looking for stdvs <= 0 indicating non-estimated units:");
##   verbose && str(verbose, value);
## 
##   verbose && exit(verbose);
## 
##   value;
## })


############################################################################
# HISTORY:
# 2007-08-17 /HB
# o Removed findUnitsTodo() from ExonChipEffectFile, because it gave the
#   same result as in the one in superclass ChipEffectFile.
# 2007-02-08 /KS
# o Created (based on SnpChipEffectFile.R following chat with HB on
#   2007-02-07).
############################################################################
