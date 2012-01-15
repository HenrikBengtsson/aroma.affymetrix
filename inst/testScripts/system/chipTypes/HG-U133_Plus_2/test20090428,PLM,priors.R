library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

csR <- AffymetrixCelSet$byName("Affymetrix-HeartBrain", 
                              chipType="HG-U133_Plus_2");

# RMA background correction
bc <- RmaBackgroundCorrection(csR);
csB <- process(bc, verbose=verbose);

# RMA quantile normalization
qn <- QuantileNormalization(csB, typesToUpdate="pm");
csN <- process(qn, verbose=verbose);

# RMA probe summarization (there are no NAs in this data set)
plm <- RmaPlm(csN);
fit(plm, verbose=verbose);

# Extract chip effects on the log2 scale
ces <- getChipEffectSet(plm);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit "another" sample based on above PLM priors
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get PLM priors
pf <- getProbeAffinityFile(plm);

# One "new" sample
csT <- extract(csN, 1);
print(csT);

# RMA probe summarization (there are no NAs in this data set)
listOfPriors <- list(probeAffinities=pf);
plmP <- RmaPlm(csT, tags="*,priors", listOfPriors=listOfPriors);
print(plmP);

# Assert that list of prior (data files) exists
priorList <- getListOfPriors(plmP);
print(priorList);

# Test to read a subset of priors
priors <- readPriorsByUnits(plmP, units=101:105);
str(priors);

fit(plmP, verbose=verbose);
cesP <- getChipEffectSet(plmP);


theta <- extractTheta(extract(ces,1), drop=TRUE);
thetaP <- extractTheta(cesP, drop=TRUE);

plot(theta, thetaP);
abline(a=0, b=1);


###########################################################################
# HISTORY:
# 2012-01-14 [HB]
# o Now sctipt also fits PLM with prior parameters.
# o Now the script is only one sample for the PLM prior part.
# 2009-04-28 [HB]
# o Added to test priors.
###########################################################################
