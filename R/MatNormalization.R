###########################################################################/**
# @RdocClass MatNormalization
#
# @title "The MatNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for systematic
#  effects in the probe intensities due to differences in the number of
#  A, C, G, and T:s and the match scores according to MAT [1].
# }
# 
# @synopsis 
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of 
#     @see "AbstractProbeSequenceNormalization".}
#   \item{unitsToFit}{The units from which the normalization curve should
#     be estimated.  If @NULL, all are considered.}
#   \item{model}{A @character string specifying the model used to fit 
#     the base-count effects.}
#   \item{numChunks}{The number of chunks to split the data into to 
#     fit the model}
#   \item{numBins}{The number of bins to use for the variance smoothing step}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \section{Requirements}{
#   This class requires that an aroma probe sequence file and aroma
#   match scores file is available for the chip type.
# }
# 
# \author{
#   Mark Robinson, WEHI.
# }
#
# \references{
#   [1] Johnson WE, Li W, Meyer CA, Gottardo R, Carroll JS, Brown M, Liu XS.
#     \emph{Model-based analysis of tiling-arrays for ChIP-chip}, PNAS, 2006.
# }
#*/###########################################################################
setConstructorS3("MatNormalization", function(..., unitsToFit=NULL, model=c("lm"), numChunks=25, numBins=200) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  extend(AbstractProbeSequenceNormalization(..., unitsToFit=unitsToFit), "MatNormalization",
    .model = model,
    .scaleResiduals = TRUE,
    .numChunks = as.integer(numChunks),
    .numBins = as.integer(numBins)
  )
})


setMethodS3("getAromaCellMatchScoreFile", "MatNormalization", function(this, ..., force=FALSE) {
  apm <- this$.apm;

  if (force || is.null(apm)) {
    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    chipType <- getChipType(cdf, fullname=FALSE);
    nbrOfCells <- nbrOfCells(cdf);
    apm <- AromaCellMatchScoreFile$byChipType(chipType, nbrOfCells=nbrOfCells, ...);
    this$.apm <- apm;
  }

  apm;
}, protected=TRUE)



setMethodS3("getAsteriskTags", "MatNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", this, collapse=collapse, ...);

  # Add model tag?
  model <- this$.model;
  tags <- c(tags, model);

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, private=TRUE)


setMethodS3("getParameters", "MatNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  params <- c(params, list(
    model = this$.model,
    numChunks = this$.numChunks
  ));

  params;
}, private=TRUE)



setMethodS3("getDesignMatrix", "MatNormalization", function(this, cells=NULL, model=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  if (is.null(cells)) {
  } else {
    # Validated below...
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving design matrix");
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchScoreFile holding match scores
  apm <- getAromaCellMatchScoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Read sequence matrix");
  sm <- readSequenceMatrix(aps, cells=cells, verbose=verbose);
  verbose && exit(verbose);
  verbose && enter(verbose, "Read match scores");
  ms <- readColumns(apm, rows=cells, verbose=verbose);
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Construct design matrix");
  nT <- rowSums(sm == "T");
  G <- (sm == "G")+0;
  A <- (sm == "A")+0;
  C <- (sm == "C")+0;
  designMatrix <- cbind(nT, A, C, G, rowSums(A)^2, rowSums(C)^2, rowSums(G)^2, nT^2, log(as.integer(ms[,1])));
  
  # Garbage collect
  rm(nT,G,A,C,ms,sm);
  gc <- gc();
  verbose && print(verbose, gc);
  
  verbose && str(verbose, designMatrix);
  verbose && cat(verbose, "object.size(designMatrix): ", 
                                             object.size(designMatrix));
  verbose && exit(verbose);

  verbose && exit(verbose);

  designMatrix;
}, private=TRUE)




setMethodS3("fitOne", "MatNormalization", function(this, df, ..., verbose=FALSE) {

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchScoreFile holding match scores
  apm <- getAromaCellMatchScoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading 'non-missing' cells to fit");
  cellsToFit <- whichVector( !(isMissing(aps) | isMissing(apm)) );
  verbose && cat(verbose, "Cells to fit:");
  verbose && str(verbose, cellsToFit);
  verbose && exit(verbose);

  numChunks <- this$.numChunks;

  # this code adopted from Dave Fourniers 17/08/2007 post to r-help mailing list
  # entitled "[R] Linear models over large datasets"

  cellsPerChunk <- ceiling(length(cellsToFit)/numChunks) + 1;
  
  verbose && enter(verbose, "Reading signals to fit");
  y <- extractMatrix(df, cells=cellsToFit, verbose=less(verbose, 10));
  verbose && exit(verbose);

  verbose && enter(verbose, "Log2 transforming signals");
  y <- log2(y);
  verbose && cat(verbose, "Target log2 probe signals:");
  verbose && str(verbose, y);
  verbose && exit(verbose);
 
  start <- xtx <- xty <- 0;
   
  while(start < length(cellsToFit)) {
  
    verbose && enter(verbose, "Working on indices over range");
    from <- start+1;
    to <- min(start+cellsPerChunk, length(cellsToFit));
    indSubset <- (from:to);
    rng <- c(from, to);
    rng <- rng / length(cellsToFit);
    verbose && cat(verbose, sprintf("[%g,%g]", rng[1], rng[2]));
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading design matrix for this iteration");
    X <- getDesignMatrix(this, cells=cellsToFit[indSubset], verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating cross products");
    xtx <- xtx + crossprod(X);
    xty <- xty + crossprod(X, y[indSubset]);
    verbose && exit(verbose);

    start <- start + cellsPerChunk;
  }
  
  verbose && enter(verbose, "Solving normal equations");
  #fit <- list(xtx=xtx, xty=xty) #,beta=solve(xtx, xty))
  verbose && exit(verbose);
  fit <- list(beta=solve(xtx, xty), scaleResiduals=this$.scaleResiduals);

  fit;
}, protected=TRUE)


setMethodS3("predictOne", "MatNormalization", function(this, fit, ..., verbose=FALSE) {

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchScoreFile holding match scores
  apm <- getAromaCellMatchScoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Allocating mu vector");
  mu <- double(nbrOfCells(aps));
  verbose && str(verbose, mu);
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Reading 'non-missing' cells to predict");
  cellsToPredict <- whichVector( !(isMissing(aps) | isMissing(apm)) );
  verbose && cat(verbose, "Cells to predict:");
  verbose && str(verbose, cellsToPredict);
  verbose && exit(verbose);

  numChunks <- this$.numChunks;

  cellsPerChunk <- ceiling(length(cellsToPredict)/numChunks) + 1;
  start <- 0;
     
  while(start < length(cellsToPredict)) {
  
    verbose && enter(verbose, "Working on indices over range");
    from <- start+1;
    to <- min(start+cellsPerChunk, length(cellsToPredict));
    indSubset <- (from:to);
    rng <- c(from, to);
    rng <- rng / length(cellsToPredict);
    verbose && cat(verbose, sprintf("[%g,%g]", rng[1], rng[2]));
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading design matrix for this iteration");
    X <- getDesignMatrix(this, cells=cellsToPredict[indSubset], verbose=verbose)
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating mu");
    mu[ cellsToPredict[indSubset] ] <- X %*% fit$beta;
    verbose && exit(verbose);

    start <- start + cellsPerChunk;
  }
  
  rm(X,indSubset);
  gc <- gc();
  
  #y <-
  
  #numBins <- this$.numBins
  #q <- quantile(mu[cellsToPredict],prob=(0:numBins)/numBins)
  #cuts<-cut(mu[cellsToPredict],breaks=q,labels=1:(length(q)-1))  # define
  #ss<-split(data.frame(resid),cuts)
  #ssvar<-sapply(ss,var)
  #v<-ssvar[as.character(cuts)]
  #for(j in 1:length(b))
  #rr<-resid/sqrt(v)

  # Return results
  mu;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod process
#
# @title "Normalizes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already normalized is re-normalized, 
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "MatNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Normalization data set for probe-sequence effects");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get (and create) the output path
  outputPath <- getPath(this);

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchScoreFile holding match scores
  apm <- getAromaCellMatchScoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading 'non-missing' cells to fit");
  cellsToFit <- whichVector( !(isMissing(aps) | isMissing(apm)) );
  verbose && cat(verbose, "Cells to fit:");
  verbose && str(verbose, cellsToFit);
  verbose && exit(verbose);

  numChunks <- this$.numChunks;
  numBins <- this$.numBins;

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize all arrays simultaneously
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- nbrOfArrays(ds);
  verbose && enter(verbose, "Normalizing ", nbrOfArrays, " arrays");
  verbose && enter(verbose, "Path: ", outputPath);

  df <- getFile(ds, 1);
  nbrOfCells <- nbrOfCells(df);
  
  xtx <- 0;
  xtyList <- as.list(double(nbrOfArrays));

  nbrOfCells <- length(cellsToFit);
  cellsPerChunk <- ceiling(nbrOfCells/numChunks) + 1;
  nbrOfChunks <- ceiling(nbrOfCells / cellsPerChunk);
  verbose && cat(verbose, "Number cells per chunk: ", cellsPerChunk);

  idxs <- 1:nbrOfCells; 
  head <- 1:cellsPerChunk;
  count <- 1;
  while(length(idxs) > 0) {
    verbose && enter(verbose, "Fitting chunk #", count, " of ", nbrOfChunks); 
    if (length(idxs) < cellsPerChunk) {
      head <- 1:length(idxs);
    }
    cc <- idxs[head];

    verbose && cat(verbose, "Cells: ");
    verbose && str(verbose, cellsToFit[cc]);
 
    verbose && enter(verbose, "Reading design matrix");
    X <- getDesignMatrix(this, cells=cellsToFit[cc], verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating cross product X'X");
    xtx <- xtx + crossprod(X);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating cross product X'y for each array");
    for (kk in seq_len(nbrOfArrays)) {
      df <- getFile(ds, kk);
      verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                          kk, getName(df), nbrOfArrays));

      y <- extractMatrix(df, cells=cellsToFit[cc], verbose=verbose);
      y <- log2(y);
      verbose && cat(verbose, "Target log2 probe signals:");
      verbose && str(verbose, y);

      xtyList[[kk]] <- xtyList[[kk]] + crossprod(X, y);

      rm(y);
      verbose && exit(verbose);
    } # for (kk ...)
    verbose && exit(verbose);

    rm(X);

    # Next chunk
    idxs <- idxs[-head]; 
    count <- count + 1;

    verbose && exit(verbose);
  }  # while (...)

  verbose && enter(verbose, "Solving for each array");
  fits <- lapply(xtyList, FUN=function(xty) {
    list(beta=solve(xtx, xty));
  });
  verbose && exit(verbose);

  # Not needed anmore
  rm(xtx, xtyList);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save model fits
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq_len(nbrOfArrays)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                              kk, getName(df), nbrOfArrays));

    fullname <- getFullName(df);

    # Store model fit 
    verbose && enter(verbose, "Saving model fit");
    # Store fit and parameters (in case someone are interested in looking
    # at them later; no promises of backward compatibility though).
    filename <- sprintf("%s,fit.RData", fullname);
    fitPathname <- Arguments$getWritablePathname(filename, 
                                                    path=outputPath, ...);
    saveObject(fits[[kk]], file=fitPathname);
    verbose && str(verbose, fits[[kk]], level=-50);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  } # for (kk ...)



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate residuals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  start <- xtx <- 0;
  mu <- vector("list", nbrOfArrays);
  cellsPerChunk <- ceiling(nbrOfCells/numChunks) + 1;
  
  while(start < nbrOfCells) {
    verbose && enter(verbose, "Working on indices over range");
    from <- start+1;
    to <- min(start+cellsPerChunk, nbrOfCells);
    indSubset <- (from:to);
    rng <- c(from, to);
    rng <- rng / nbrOfCells;
    verbose && cat(verbose, sprintf("[%g,%g]", rng[1], rng[2]));
    verbose && exit(verbose);

    verbose && enter(verbose, "Set of cells");
    cellsChunk <- cellsToFit[indSubset];
    verbose && str(verbose, cellsChunk);
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading design matrix for this iteration");
    X <- getDesignMatrix(this, cells=cellsChunk, verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating model fits");
    xtx <- xtx + crossprod(X);

    verbose && enter(verbose, "Processing ", nbrOfArrays, " arrays");
    for (kk in seq_len(nbrOfArrays)) {
      df <- getFile(ds, kk);
      verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                              kk, getName(df), nbrOfArrays));

      mu <- X %*% fits[[kk]]$beta;
      mu <- as.double(mu);
      verbose && str(verbose, mu);

      fullname <- getFullName(df);
      filename <- sprintf("%s.CEL", fullname);
      pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");
      createFrom(df, filename=pathname, path=NULL, verbose=less(verbose));
      verbose && exit(verbose);

      updateCel(pathname, indices=cellsChunk, intensities=2^mu, verbose=TRUE);
      verbose && exit(verbose);
    } # for (kk ...)
    verbose && exit(verbose);

    rm(mu);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);

    start <- start + cellsPerChunk;
    #start <- nbrOfCells + 1;

  } # while (start < ...)



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Scale residuals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq_len(nbrOfArrays)) {
      verbose && enter(verbose, "Binning predicted values, calculating and scaling residuals");
      df <- getFile(ds, kk);

      y <- extractMatrix(df, cells=cellsToFit, verbose=verbose);
      y <- log2(y);

      fullname <- getFullName(df);
      filename <- sprintf("%s.CEL", fullname);
      pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

      mu <- readCel(pathname, indices=cellsToFit, readOutliers=FALSE, readHeader=FALSE, readMasked=FALSE, verbose=less(verbose,10))$intensities;
      mu <- log2(mu);
      r <- y - mu;

      q <- quantile(mu, probs=(0:numBins)/numBins);
      cuts <- cut(mu, breaks=q, labels=1:(length(q)-1));  # define
      ss <- split(r, cuts);
      ssvar <- sapply(ss, FUN=var);
      v <- ssvar[as.character(cuts)];
      r <- r / sqrt(v);
      r <- as.double(r);

      #return(list(y=y,mu=mu,r=r))

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");
      createFrom(df, filename=pathname, path=NULL, verbose=less(verbose));
      verbose && exit(verbose);

      updateCel(pathname, indices=cellsToFit, intensities=2^r, verbose=TRUE);
      rm(q,ss,ssvar,v,r,y);
      gc <- gc();
      verbose && print(verbose, gc);

      verbose && exit(verbose);
  } # for (kk ...)

  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);
  
  invisible(outputDataSet);
})



############################################################################
# HISTORY:
# 2009-05-23 [HB]
# o Updated some minor format mistakes in the verbose output.
# 2008-11-28 [HB]
# o Added protected getCrossProductXTX().  Still not used.
# o Modified the first loop over chunks that calculated cross products such
#   that it is constant in number of arrays.  The processing over chunks is
#   now also done as we do it elsewhere in the package.
# o fitOne() and predictOne() are never used?!?
# o Updated the Rdocs.
# o Renamed getAromaCellMatchscoreFile() to getAromaCellMatchScoreFile().
# 2008-10-29 [MR]
# o Created from BaseCountNormalization.R
############################################################################
