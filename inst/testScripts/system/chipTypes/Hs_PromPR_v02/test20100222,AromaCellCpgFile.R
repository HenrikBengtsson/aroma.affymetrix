###########################################################################/**
# File type test
#
# Description:
# This test verifies that aroma.affymetrix can create AromaCellCpgFile
# and AromaCellPositionFile objects for the (promoter) tiling array.
#
# Author: Mark Robinson
# Created: 2010-02-22
# Last modified: 2010-03-14
#
#*/###########################################################################

library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup the chip type
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02");
print(cdf);

cdfU <- getUniqueCdf(cdf);
print(cdfU);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test allocation, writing and reading of 'acp' object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acp <- AromaCellPositionFile$allocateFromCdf(cdfU, tags=c("*", "test"), 
                                            overwrite=TRUE, verbose=verbose);
print(acp);

nRandCells <- 20;

cells <- sample(seq_len(nbrOfCells(cdfU)), nRandCells);

chRand <- sample(seq_len(22), nRandCells);
posRand <- sample(seq_len(1e6), nRandCells);

acp[cells,1] <- chRand;
acp[cells,2] <- posRand;

tmpMatrix <- acp[cells,];

stopifnot((tmpMatrix[,1] == chRand) & (tmpMatrix[,2] == posRand));

rm(chRand, posRand, tmpMatrix);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test allocation, writing and reading of 'acc' object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AromaCellCpgFile$allocateFromCdf(cdfU, tags=c("*", "test"), 
                                            overwrite=TRUE, verbose=verbose);
print(acc);

cells <- sample(seq_len(nbrOfCells(cdfU)), nRandCells);
cpgRand <- rnorm(nRandCells);

acc[cells,1] <- 2^cpgRand;

tmpMatrix <- acc[cells,];

ss <- sum( (log2(tmpMatrix[,1])-cpgRand)^2 );

stopifnot(ss < 1e-8);
