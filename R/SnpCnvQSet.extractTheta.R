setMethodS3("extractTheta", "SnpCnvQSet", function(this, ..., transform=function(y, ...) { 2^y }, addNames=TRUE, verbose=FALSE) {
  eSet <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  nbrOfUnits <- nrow(eSet);
  nbrOfSamples <- ncol(eSet);
  nbrOfGroups <- 2;  # (thetaA, thetaB)

  # Extract sample names
  sampleNames <- sampleNames(eSet);
  sampleNames <- gsub("[.](cel|CEL)$", "", sampleNames);
  sampleNames <- gsub(",.*$", "", sampleNames);

  # Extract unit names
  snpNames <- featureNames(eSet);

  # Allocate result object
  naValue <- as.double(NA);
  theta <- array(naValue, dim=c(nbrOfUnits, nbrOfGroups, nbrOfSamples));
  dimnames(theta)[[3]] <- sampleNames;
  if (addNames)
    dimnames(theta)[[1]] <- snpNames;

  # Populate with estimates
  theta[,1,] <- transform(oligo:::thetaA(eSet));
  theta[,2,] <- transform(oligo:::thetaB(eSet));

  theta;
})


############################################################################
# HISTORY:
# 2008-12-05
# o Created.
############################################################################
