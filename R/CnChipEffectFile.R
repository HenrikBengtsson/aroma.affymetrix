###########################################################################/**
# @RdocClass CnChipEffectFile
#
# @title "The CnChipEffectFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in a copy-number probe-level
#  models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "SnpChipEffectFile".}
#   \item{combineAlleles}{A @logical indicating if the signals from allele A 
#     and allele B are combined or not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#   An object of this class is typically part of a @see "CnChipEffectSet".
# }
#*/###########################################################################
setConstructorS3("CnChipEffectFile", function(..., combineAlleles=FALSE) {
  this <- extend(SnpChipEffectFile(...), "CnChipEffectFile",
    combineAlleles = combineAlleles
  );

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})

setMethodS3("getParameters", "CnChipEffectFile", function(this, ...) {
  params <- NextMethod(generic="getParameters", object=this, ...);
  params$combineAlleles <- this$combineAlleles;
  params;
})


setMethodS3("getCellIndices", "CnChipEffectFile", function(this, units=NULL, ..., force=FALSE, .cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  # Argument 'units':
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "getCellIndices.CnChipEffectFile()");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force || .cache) {
    chipType <- getChipType(getCdf(this));
    params <- getParameters(this);
    key <- list(method="getCellIndices", class=class(this)[1], 
                chipType=chipType, params=params, units=units, ...);
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
  # Get units in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(units))
    units <- seq(length=nbrOfUnits(cdf));

  cells <- lapplyInChunks(units, function(unitChunk) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get and restructure cell indices
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## NOTE: NextMethod() does not work from within another function
##    cells <- NextMethod("getCellIndices", this, units=unitChunk, ..., 
##                            force=force, .cache=FALSE, verbose=verbose);
    cells <- getCellIndices.SnpChipEffectFile(this, units=unitChunk, ..., 
                             force=force, .cache=FALSE, verbose=verbose);
    gc <- gc();
    verbose && print(verbose, gc);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Combine alleles?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # If combining alleles, return only every second group.
    # In order to improve readability we merge the names of alleles groups
    # combined, e.g. groups 'C' and 'G' become group 'CG'.
    if (this$combineAlleles) {
      verbose && enter(verbose, "Combining alleles");
      # Hard-wiring 1, 2 & 4 groups speed things up 3 times!

      cells <- applyCdfGroups(cells, function(groups) {
        ngroups <- length(groups);
        names <- names(groups);
        if (ngroups == 4) {
          groups <- .subset(groups, c(1,3));
          names <- paste(.subset(names, c(1,3)), .subset(names, c(2,4)), sep="");
        } else if (ngroups == 2) {
          groups <- .subset(groups, 1);
          names <- paste(.subset(names, 1), .subset(names, 2), sep="");
        } else if (ngroups == 1) {
          groups <- .subset(groups, 1);
        } else {
          odds <- seq(from=1, to=ngroups, by=2);
          evens <- seq(from=2, to=ngroups, by=2);
          groups <- .subset(groups, odds);
          names <- paste(.subset(names, odds), .subset(names, evens), sep="");
        }
        names(groups) <- names;
        groups;
      }) # applyCdfGroups()
      verbose && printf(verbose, "Number of units: %d\n", length(cells));
      verbose && exit(verbose);
    }

    cells;
  }, chunkSize=100e3, verbose=less(verbose)) # lapplyInChunks()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.cache) {
    # In-memory or on-file cache?
    size <- object.size(cells);
      verbose && printf(verbose, "Object size: %.1fMB\n", size/1024^2);
    if (size < 10e6) { 
      # In-memory cache for objects < 10Mb.
      this$.cellIndices <- list();
      this$.cellIndices[[id]] <- cells;
      verbose && cat(verbose, "Result cached in memory");
    } else {
      # On-file cache
      # Keep, in-memory cache.
      if (!is.list(this$.cellIndices))
        this$.cellIndices <- list();
      this$.cellIndices[[id]] <- NULL;
      saveCache(cells, key=list(id), dirs=dirs);
      verbose && cat(verbose, "Result cached to file");
    }
  }

  verbose && exit(verbose);

  cells;
})


setMethodS3("readUnits", "CnChipEffectFile", function(this, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Check for cached data
  key <- list(method="readUnits", class=class(this)[1],
              pathname=getPathname(this),
              combineAlleles=this$combineAlleles, ...);
  id <- digest2(key);
  res <- this$.readUnitsCache[[id]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "readUnits.CnChipEffectFile(): Returning cached data");
    return(res);
  }

  # Retrieve the data
  res <- NextMethod("readUnits", this, ..., force=TRUE, cache=FALSE, verbose=verbose);


  # Store read units in cache?
  if (cache) {
    verbose && cat(verbose, "readUnits.CnChipEffectFile(): Updating cache");
    this$.readUnitsCache <- list();
    this$.readUnitsCache[[id]] <- res;
  }

  res;
})


setMethodS3("mergeStrands", "CnChipEffectFile", function(this, ...) {
  cfM <- NextMethod("mergeStrands", this, ...);
  cfM$combineAlleles <- this$combineAlleles;
  cfM;
})

############################################################################
# HISTORY:
# 2008-02-19
# o BUG FIX: getCellIndices() of CnChipEffectFile would return an error
#   if 'units==NULL'.
# 2007-03-01
# o BUG FIX: getCellIndices() would give "Error in fcn(.subset2(unit, 
#   "groups"), ...) : object "odds" not found" for units with other than
#   1, 2, or 4 groups.
# 2007-01-20
# o Added mergeStrands().
# 2007-01-06
# o Now getCellIndices() caches large objects to file and small in memory.
# o Made getCellIndices() three times faster by some hardwired code.
# 2006-11-28
# o Added readUnits() to override caching mechanism of superclasses.
# 2006-09-20
# o BUG FIX: Typo. Remove an argument but tried to use inside.
# 2006-09-17
# o Added an in-memory cache for getCellIndices().
# 2006-09-12
# o Updated and probably working. When combining alleles, the names of the
#   groups returned consist of the allele A and allele group names.
# 2006-09-11
# o Created.
############################################################################
