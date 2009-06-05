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
# With priors
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
pf <- getProbeAffinityFile(plm);

# RMA probe summarization (there are no NAs in this data set)
plmP <- RmaPlm(csN);
setListOfPriors(plmP, list(probeAffinities=pf));

# Assert that list of prior (data files) exists
priorList <- getListOfPriors(plmP);
print(priorList);

# Test to read a subset of priors
priors <- readPriorsByUnits(plmP, units=101:105);
str(priors);


###########################################################################
# HISTORY:
# 2009-04-28 [HB]
# o Added to test priors.
###########################################################################
