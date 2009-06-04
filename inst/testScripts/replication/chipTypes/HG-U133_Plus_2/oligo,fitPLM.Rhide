###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the RMA 
# chip-effect estimates as estimated by oligo.
# The setup is the same as in affyPLM,fitPLM.R.
#
# Author: Henrik Bengtsson
# Created: 2008-12-04 (from affyPLM,fitPLM.R)
# Last modified: 2008-12-04
###########################################################################

library("aroma.affymetrix");
library("oligo");

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

doPlot <- TRUE;
saveImg <- TRUE;

# ----------------------------------
# RMA estimates by aroma.affymetrix 
# ----------------------------------
verbose && enter(verbose, "RMA by aroma.affymetrix");

csR <- AffymetrixCelSet$byName("Affymetrix-HeartBrain", 
                              chipType="HG-U133_Plus_2");

# RMA background correction
bc <- RmaBackgroundCorrection(csR);
csB <- process(bc, verbose=verbose);

# RMA quantile normalization
qn <- QuantileNormalization(csB, typesToUpdate="pm");
csN <- process(qn, verbose=verbose);

# RMA probe summarization (there are no NAs in this data set)
plm <- RmaPlm(csN, flavor="oligo");
fit(plm, verbose=verbose);

# Extract chip effects on the log2 scale
ces <- getChipEffectSet(plm);
theta <- extractMatrix(ces, returnUgcMap=TRUE);
theta <- log2(theta);
ugcMap <- attr(theta, "unitGroupCellMap");
rownames(theta) <- getUnitNames(getCdf(ces), ugcMap[,"unit"]);

verbose && exit(verbose);


# ----------------------
# RMA estimates by oligo
# ----------------------
verbose && enter(verbose, "RMA by oligo");
verbose && print(verbose, sessionInfo());

library("pd.hg.u133.plus.2");

raw <- read.celfiles(filenames=getPathnames(csR));
eSet <- rma(raw);
theta0 <- exprs(eSet);

verbose && exit(verbose);


# --------------------------------
# Compare the two implementations
# --------------------------------
# Reorder the aroma.affymetrix estimates
o <- match(rownames(theta0), rownames(theta));
theta <- theta[o,];

# (a) Assert correlations
cors <- sapply(1:ncol(theta), FUN=function(cc) cor(theta[,cc], theta0[,cc]));
print(cors);
print(range(cors));
stopifnot(all(cors > 0.99995));

# (b) Assert differences
e <- (theta - theta0);
stopifnot(mean(as.vector(e^2)) < 1e-3);
stopifnot(sd(as.vector(e^2)) < 1e-3);
stopifnot(quantile(abs(e), 0.99) < 0.05);
stopifnot(max(abs(e)) < 0.085);

if (doPlot) {
  if (saveImg) {
    pngDev <- findPngDevice();
    devNew("pngDev", "replication-oligo,fitPLM.png", width=640, height=640);
  }

  layout(matrix(1:9, ncol=3, byrow=TRUE));

  xlab <- expression(log[2](theta[oligo]));
  ylab <- expression(log[2](theta[aroma.affymetrix]));
  for (kk in seq(length=ncol(theta))) {
    main <- colnames(theta)[kk];
    plot(theta0[,kk], theta[,kk], pch=".", xlab=xlab, ylab=ylab, main=main);
    abline(0,1, col="blue");
  }

  xlab <- expression(log[2](theta[aroma.affymetrix]/theta[oligo]));
  plotDensity(e, xlab=xlab);

  devDone();
}

###########################################################################
# HISTORY:
# 2008-12-04 [HB]
# o Created, but not tested because I miss package 'pd.hg.u133.plus.2'.
###########################################################################
