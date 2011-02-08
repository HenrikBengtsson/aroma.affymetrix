library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, verbose=log);
keep <- 1:6;
csR <- extract(csR, keep);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR);
print(acc);
csC <- process(acc, verbose=log);
print(csC);
stopifnot(identical(getNames(csC), getNames(csR)));



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for allele-specific CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=FALSE, shift=300);
print(plm);

fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);
stopifnot(identical(getNames(ces), getNames(csR)));



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plotting allele B frequences and raw copy numbers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- getCdf(ces);
gi <- getGenomeInformation(cdf);
units <- getUnitsOnChromosome(gi, 2);
data <- extractDataFrame(ces, units=units, addNames=TRUE, verbose=log);
data[,"x"] <- getPositions(gi, data[,"unit"]);

keep <- match(getNames(ces), colnames(data));
theta <- as.matrix(data[,keep]);
data <- data[,-keep];

isA <- (data$group==1);
isB <- (data$group==2);
thetaA <- theta[isA,];
thetaB <- theta[isB,];

# Total raw copy numbers
theta <- thetaA + thetaB;
thetaR <- rowMedians(theta, na.rm=TRUE);
M <- log2(theta/thetaR);

# Allele B frequences
B <- thetaB/theta;

# Position (in Mb)
x <- data[isA,"x"] / 1e6;

layout(matrix(1:2, nrow=2));
par(mar=c(5,4,2,2)+0.1);
xlim <- range(x, na.rm=TRUE);
xlab <- "Physical position (in Mb)";
Blab <- expression(theta[B]/theta);
Mlab <- expression(log[2](theta/theta[R]));

cc <- 1;
plot(NA, xlim=xlim, ylim=c(0,1), xlab=xlab, ylab=Blab);
abline(h=1/2, col="#cccccc");
points(x,B[,cc], pch=".", cex=2);
plot(NA, xlim=xlim, ylim=c(-1,1)*3, xlab=xlab, ylab=Mlab);
abline(h=0, col="#cccccc");
points(x,M[,cc], pch=".", cex=2);

devDone();

