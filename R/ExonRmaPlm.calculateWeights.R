
setMethodS3("calculateWeights", "ExonRmaPlm", function(this, units=NULL, ram=NULL, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  resFcn <- function(unit, mergeGroups) {
    nbrOfGroups <- length(unit);
    if (mergeGroups) {
      y <- do.call("rbind", base::lapply(unit, .subset2, "eps"));
      y <- log2(y);
      madMerged <- 1.4826 * median(abs(y));
    }
    res <- base::lapply(1:nbrOfGroups, FUN=function(gg) {
      y <- .subset2(.subset2(unit, gg), "eps");
      y <- log2(y);
      mad <- 1.4826 * median(abs(y));
      if (mad==0) {
        return(matrix(data=1, nrow=nrow(y), ncol=ncol(y)));        
      }
      if (mergeGroups) {
        return(matrix(MASS::psi.huber(y/madMerged), ncol=ncol(y)));
      } else {
        return(matrix(MASS::psi.huber(y/mad), ncol=ncol(y)));
      }
    })
    res;
  } # resFcn()
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  rs <- calculateResiduals(this, verbose=verbose);
  ws <- getWeightsSet(this, verbose=verbose);
  nbrOfArrays <- nbrOfArrays(rs);

  verbose && enter(verbose, "Calculating PLM weights");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get data and parameter objects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this);
  if (is.null(ds)) {
    throw("No data set specified for PLM: ", getFullName(this));
  }

  cdf <- getCdf(ds);
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
  } else {
    nbrOfUnits <- length(units);
  }

  if (force) {
    unitsToDo <- units;
  } else {
    unitsToDo <- findUnitsTodo(ws, units=units);
  }
  
  verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits);

  unitsPerChunk <- ram * 100000/nbrOfArrays(getDataSet(this));
  unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range=c(1,Inf));
  
  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && printf(verbose, "Number of chunks: %d (%d units/chunk)\n",
                    nbrOfChunks, unitsPerChunk);
  head <- 1:unitsPerChunk;

  verbose && enter(verbose, "Extracting unit data");
  count <- 1;
  while (length(unitsToDo) > 0) {
    if (length(unitsToDo) < unitsPerChunk) {
      head <- 1:length(unitsToDo);
    }
    units <- unitsToDo[head];
    verbose && printf(verbose, "Chunk #%d of %d (%d units)\n",
                                        count, nbrOfChunks, length(units));

    residualsList <- readUnits(rs, units=units, verbose=less(verbose), stratifyBy="pm");

    verbose && enter(verbose, "Calculating weights");
    weightsList <- base::lapply(residualsList, FUN=resFcn, mergeGroups=this$mergeGroups);
    verbose && exit(verbose);
    
# update output files
    
    verbose && enter(verbose, "Storing weights");

    cdf <- getCellIndices(getCdf(ds), units=units, stratifyBy="pm", ...);
    
    for (kk in seq(ds)) {
      wf <- getFile(ws, kk);

      verbose && enter(verbose, sprintf("Array #%d ('%s')", kk, getName(wf)));

      data <- base::lapply(weightsList, function(unit) {
        base::lapply(unit, function(group) {
          nrow <- nrow(group); 
          list(
               intensities=2^group[,kk], 
               stdvs=rep(1, nrow), 
               pixels=rep(1, nrow)
               );
        })
      });
      
      updateCelUnits(getPathname(wf), cdf=cdf, data=data);

      verbose && exit(verbose);
    }
    
    verbose && exit(verbose);

    unitsToDo <- unitsToDo[-head];
    count <- count + 1;
  }

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);
  
  verbose && exit(verbose);

  invisible(ws);
})


##########################################################################
# HISTORY:
# 2007-02-15 
# o Based on ProbeLevelModel.calculateResiduals
#   and QualityAssessmentModel.getWeights
##########################################################################
