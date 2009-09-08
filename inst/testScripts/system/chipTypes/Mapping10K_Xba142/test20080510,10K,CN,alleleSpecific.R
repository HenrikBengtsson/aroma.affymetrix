library("aroma.affymetrix");
library("geneplotter");

log <- Arguments$getVerbose(-4, timestamp=TRUE);



pngDev <- findPngDevice();
imgFormat <- "png";
figPath <- "figures";


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
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaSnpPlm(csC, mergeStrands=TRUE, shift=300);
print(plm);

fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);
stopifnot(identical(getNames(ces), getNames(cs)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extraction test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract unit and group name and indices together with chip effects
data <- extractDataFrame(ces, units=1:50, addNames=TRUE, verbose=log);
print(data[1:50,1:8]);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract (thetaA, thetaB)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
theta <- extractTheta(ces, groups=1:2, verbose=log);
dimnames(theta)[[2]] <- c("A", "B");
thetaA <- theta[,"A",];
thetaB <- theta[,"B",];
theta <- thetaA + thetaB;
freqB <- thetaB / theta;


tlim <- c(8,15);
Blim <- c(0,1);
tlab <- expression(log[2](theta));
Blab <- expression(theta[B]/theta);

nbrOfArrays <- ncol(theta);
layout(matrix(1:nbrOfArrays, ncol=2, byrow=TRUE));
par(mar=c(3.8,4,3,1)+0.1);
for (kk in seq(length=nbrOfArrays)) {
  name <- colnames(theta)[kk];
  smoothScatter(log2(theta[,kk]), freqB[,kk], 
                xlim=tlim, ylim=Blim, xlab=tlab, ylab=Blab, main=name);
}
devDone();


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract (log2(total), freqB)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- extractTotalAndFreqB(ces, verbose=log);
data[,1,] <- log2(data[,1,]);

tlim <- c(9,15);

upperPanel <- function(x,y, pch=".", ...) {
  xy <- data[,1,c(x,y)];
  xy <- (xy - tlim[1])/diff(tlim);
  abline(a=0,b=1, col="#999999");
  points(xy, pch=pch);
}

lowerPanel <- function(x,y, pch=".", ...) {
  xy <- data[,2,c(x,y)];
  abline(a=0,b=1, col="#999999");
  abline(a=1,b=-1, col="#999999");
  points(xy, pch=pch);
}

diagPanel <- function(x, pch=".", ...) {
  xy <- data[,,x];
  xy[,1] <- (xy[,1] - tlim[1])/diff(tlim);
  smoothScatter(xy, add=TRUE);
  abline(h=1/2, col="#999999");
}

pairs <- matrix(1:nbrOfArrays, nrow=1, ncol=nbrOfArrays);
colnames(pairs) <- dimnames(data)[[length(dim)]];

if (imgFormat == "png") {
  filename <- sprintf("%s,pairs.png", getFullName(ces));
  pathname <- filePath(figPath, filename);
  devNew("pngDev", pathname, width=1024, height=1024);
}
pairs(pairs, upper.panel=upperPanel, lower.panel=lowerPanel, 
             diag.panel=diagPanel, xlim=c(0,1), ylim=c(0,1));
devDone();
