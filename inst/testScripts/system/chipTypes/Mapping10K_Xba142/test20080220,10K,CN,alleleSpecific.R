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

# Extract (thetaA, thetaB)
cdf <- getCdf(ces);
snps <- indexOf(cdf, "^SNP");
theta <- extractMatrix(ces, units=snps, returnUgcMap=TRUE, verbose=log);
ugcMap <- attr(theta, "unitGroupCellMap");
theta <- c(theta[ugcMap$group==1,], theta[ugcMap$group==2,]);
dimnames <- list(getUnitNames(cdf, units=snps), c("A", "B"), getNames(ces));
dim <- sapply(dimnames, FUN=length);
theta <- array(theta, dim=dim, dimnames=dimnames);

ltheta <- log2(theta);
a <- (ltheta[,"A",]+ltheta[,"B",])/2;
m <- ltheta[,"A",]-ltheta[,"B",];

Alim <- c(5,14); Mlim <- diff(Alim)*c(-1,1);
Alab <- expression(1/2%*%log[2](theta[A]*theta[B]));
Mlab <- expression(log[2](theta[A]/theta[B]));

layout(matrix(1:9, ncol=3, byrow=TRUE));
par(mar=c(3.8,4,3,1)+0.1);
for (kk in 1:nbrOfArrays(ces)) {
  name <- getNames(ces)[kk];
  plot(a[,kk], m[,kk], pch=".", xlim=Alim, ylim=Mlim, xlab=Alab, ylab=Mlab, main=name);
}
devDone();
