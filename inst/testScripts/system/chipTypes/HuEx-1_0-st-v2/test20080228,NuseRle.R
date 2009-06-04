library("aroma.affymetrix");

log <- Arguments$getVerbose(-3, timestamp=TRUE);



# The reason for observing a small difference is scores below,
# is that when saving to file the estimates are rounded of to
# floats, whereas the estimates calculated in memory are kept
# in full precision.
tol <- 1e-5;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set (from previous file)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
name <- "Affymetrix-HeartBrain";
chipType <- "HuEx-1_0-st-v2";
cdf <- AffymetrixCdfFile$byChipType(chipType, 
                                         tags="coreR3,A20071112,EP");
cs <- AffymetrixCelSet$byName(name=name, cdf=cdf);

# Background correction and normalization
bc <- RmaBackgroundCorrection(cs);
csBC <- process(bc, verbose=log);
qn <- QuantileNormalization(csBC, typesToUpdate="pm");
csN <- process(qn, verbose=log);

# Probe-level summarization
plmList <- list(
  merge   = ExonRmaPlm(csN, mergeGroups=TRUE),
  noMerge = ExonRmaPlm(csN, mergeGroups=FALSE)
);

cesList <- lapply(plmList, FUN=getChipEffectSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RLE/NUSE for first 100 units
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
units <- 1:100;

# Make sure the units are fitted
dummy <- lapply(plmList, FUN=fit, units=units, verbose=log);
rm(dummy);

for (name in names(cesList)) {
  log && enter(log, name);

  ces <- cesList[[name]];

  # Assert correctness of unit-specific RLE scores
  theta <- extractMatrix(ces, field="theta", units=units);
  thetaR <- 2^rowMedians(log2(theta), na.rm=TRUE);
  rle0 <- log2(theta/thetaR);
  rle1 <- extractMatrix(ces, field="RLE", units=units, verbose=log);
  stopifnot(all.equal(rle1, rle0, tolerance=tol));

  # Assert correctness of boxplot statistics of RLE scores
  stats0 <- boxplot.stats(rle0[,1]);
  stats1 <- boxplotStats(ces, type="RLE", arrays=1:2, subset=units);
  stopifnot(all.equal(stats1[[1]], stats0, tolerance=tol));
  
  # Assert correctness of unit-specific NUSE scores
  se <- extractMatrix(ces, field="sdTheta", units=units);
  seR <- 2^rowMedians(log2(se), na.rm=TRUE);
  nuse0 <- log2(se)/log2(seR);
  nuse1 <- extractMatrix(ces, field="NUSE", units=units, verbose=log);
  stopifnot(all.equal(nuse1, nuse0, tolerance=tol));

  # Assert correctness of boxplot statistics of NUSE scores
  stats0 <- boxplot.stats(nuse0[,1]);
  stats1 <- boxplotStats(ces, type="NUSE", arrays=1:2, subset=units);
  stopifnot(all.equal(stats1[[1]], stats0, tolerance=tol));

  log && exit(log);
} # for (name in ...)
