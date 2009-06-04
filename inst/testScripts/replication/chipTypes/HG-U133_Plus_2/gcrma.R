###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the gcRMA 
# chip-effect estimates as estimated by gcrma.
#
# Author: Mark Robinson and Henrik Bengtsson
# Created: 2009-05-17
# Last modified: 2009-05-17
#
# Data set:
#  rawData/
#   Affymetrix-HeartBrain/
#    HG-U133_Plus_2/
#     u1332plus_ivt_cerebellum_A.CEL [13555904 bytes]
#     u1332plus_ivt_cerebellum_B.CEL [13550687 bytes]
#     u1332plus_ivt_cerebellum_C.CEL [13551860 bytes]
#     u1332plus_ivt_heart_A.CEL      [13554731 bytes]
#     u1332plus_ivt_heart_B.CEL      [13553255 bytes]
#     u1332plus_ivt_heart_C.CEL      [13551203 bytes]
#  Source: Affymetrix Tissue samples, 2007.  http://www.affymetrix.com/
#  support/technical/sample_data/hugene_1_0_array_data.affx
###########################################################################

library("aroma.affymetrix");
library("gcrma");  # gcrma()

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

doPlot <- TRUE;
saveImg <- TRUE;

# ----------------------------------
# RMA estimates by aroma.affymetrix 
# ----------------------------------
verbose && enter(verbose, "gcRMA by aroma.affymetrix");

cdf <- AffymetrixCdfFile$byChipType("HG-U133_Plus_2");
csR <- AffymetrixCelSet$byName("Affymetrix-HeartBrain", cdf=cdf);

# RMA background correction
bc <- GcRmaBackgroundCorrection(csR);
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
rownames(theta) <- getUnitNames(cdf);
verbose && str(verbose, theta);

verbose && exit(verbose);


# ------------------------
# gcRMA estimates by gcrma
# ------------------------
verbose && enter(verbose, "gcRMA by gcrma");
verbose && print(verbose, sessionInfo());

raw <- ReadAffy(filenames=getPathnames(csR));
verbose && print(verbose, raw);

es <- gcrma(raw, verbose=TRUE);
verbose && print(verbose, es);

theta0 <- exprs(es);
verbose && str(verbose, theta0);
verbose && exit(verbose);


# ---------------------------------
# Comparing the two implementations
# ---------------------------------
verbose && enter(verbose, "Comparing the two implementations");
# Reorder the aroma.affymetrix estimates
o <- match(rownames(theta0), rownames(theta));
theta <- theta[o,];

# (a) Assert correlations
cors <- sapply(1:ncol(theta), FUN=function(cc) cor(theta[,cc], theta0[,cc]));
print(cors);
print(range(cors));
stopifnot(all(cors > 0.999));

# (b) Assert differences
e <- (theta - theta0);
stopifnot(mean(as.vector(e^2)) < 0.003);
stopifnot(sd(as.vector(e^2)) < 0.02);
stopifnot(quantile(abs(e), 0.99) < 0.20);
stopifnot(max(abs(e)) < 1.5);
verbose && exit(verbose);

if (doPlot) {
  if (saveImg) {
    pngDev <- findPngDevice();
    devNew("pngDev", "replication-gcrma.png", width=640, height=640);
  }

  layout(matrix(1:9, ncol=3, byrow=TRUE));

  xlab <- expression(log[2](theta[gcrma]));
  ylab <- expression(log[2](theta[aroma.affymetrix]));
  for (kk in seq(length=ncol(theta))) {
    main <- colnames(theta)[kk];
    plot(theta0[,kk], theta[,kk], pch=".", xlab=xlab, ylab=ylab, main=main);
    abline(0,1, col="blue");
  }

  xlab <- expression(log[2](theta[aroma.affymetrix]/theta[gcrma]));
  plotDensity(e, xlab=xlab);

  devDone();
}

###########################################################################
# HISTORY:
# 2009-05-17 [HB]
# o Adopted from Mark Robinson's test script.
###########################################################################
