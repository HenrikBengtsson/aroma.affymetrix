library("aroma.affymetrix");
log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, verbose=log);
keep <- 1:6;
csR <- extract(csR, keep);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR);
print(acc);
csC <- process(acc, verbose=log);
print(csC);
stopifnot(identical(getNames(csC), getNames(csR)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=TRUE, shift=300);
print(plm);

fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);
stopifnot(identical(getNames(ces), getNames(csR)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization (toward an average effect)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces);
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);
stopifnot(identical(getNames(cesN), getNames(ces)));

theta <- extractTheta(cesN, drop=TRUE);
thetaR <- rowMedians(theta, na.rm=TRUE);
M <- log2(theta/thetaR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization (toward a constant effect)
# Note, this a pure single-array approach.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces, target="zero", tags="*,z");
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);
stopifnot(identical(getNames(cesN), getNames(ces)));

theta <- extractTheta(cesN, drop=TRUE);
thetaR <- rowMedians(theta, na.rm=TRUE);
M2 <- log2(theta/thetaR);

# When calculating the log-ratios, the above two approaches should
# give equals results, because the effects should cancel out regardless.
stopifnot(mean(abs(M2-M), na.rm=TRUE) < 1e-3);
