library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);



dataSetName <- "Jeremy_2007-10k";
chipType <- "Mapping10K_Xba142";

# Expected sample names
sampleNames <- c("0001-7", "0002-10", "0004-13", "0005-14", "0007-18", 
                      "0008-19", "0010-22", "2-DPrrr", "MH12", "MH18");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
keep <- 1:6;
cs <- extract(cs, keep);
sampleNames <- sampleNames[keep];
print(cs);
stopifnot(identical(getNames(cs), sampleNames));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(cs);
print(acc);
csC <- process(acc, verbose=log);
print(csC);
stopifnot(identical(getNames(csC), getNames(cs)));



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for allele-specific CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=FALSE, shift=300);
print(plm);

fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);
stopifnot(identical(getNames(ces), getNames(cs)));



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

