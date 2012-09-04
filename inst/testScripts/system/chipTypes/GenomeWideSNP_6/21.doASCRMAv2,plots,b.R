library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-50, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE13372,testset";
chipType <- "GenomeWideSNP_6,Full";

csR <- AffymetrixCelSet$byName(dataSetName, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract (thetaA, thetaB) for copy-neutral chromosomes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesN <- res$cesN;
cdf <- getCdf(cesN);
gi <- getGenomeInformation(cdf);
units <- getUnitsOnChromosomes(gi, 1:22);
theta <- extractTheta(cesN, units=units);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Estimate copy numbers on the natural scale
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thetaR <- rowMedians(theta[,1,]+theta[,2,], na.rm=TRUE);
C <- 2*(theta[,1,]+theta[,2,])/thetaR;
CA <- 2*theta[,1,]/thetaR;
CB <- 2*theta[,2,]/thetaR;
freqB <- CB/C;

mu <- colMedians(C, na.rm=TRUE);
muA <- colMedians(CA, na.rm=TRUE);
muB <- colMedians(CB, na.rm=TRUE);

Clim <- c(-0.2,3.2);
xlab <- "Freq B (thetaB/theta)";
Clab <- "Copy number";
CAlab <- "Copy number (Allele A)";
CBlab <- "Copy number (Allele B)";

fig <- 1;
if (!devIsOpen(fig <- fig + 1)) {
  devNew();
  layout(matrix(1:9, nrow=3, byrow=TRUE));
  par(mar=c(4,4,0.5,0.5)+0.1);
  centers <- matrix(c(0,2, 1,1, 2,0), nrow=3, ncol=2, byrow=TRUE);
  for (cc in 1:ncol(C)) {
    xx <- CA[,cc];
    yy <- CB[,cc];
    X <- cbind(xx,yy);
    ok <- (is.finite(X) & -1 < X & X < 30);
    ok <- ok[,1] & ok[,2];
    X <- X[ok,];
    smoothScatter(X, xlim=Clim, ylim=Clim, xlab=CAlab, ylab=CBlab);
    abline(h=0:3, lty=3, lwd=1);
    abline(v=0:3, lty=3, lwd=1);
  
    # Plot centers
    fit <- kmeans(X, centers=centers);
    print(fit$centers);
    points(fit$centers, pch=19, col="red");
  }
  devDone();
}

if (!devIsOpen(fig <- fig + 1)) {
  devNew();
  layout(matrix(1:9, nrow=3, byrow=TRUE));
  par(mar=c(4,4,0.5,0.5)+0.1);
  centers <- matrix(c(0,2, 1/2,2, 1,2), nrow=3, ncol=2, byrow=TRUE);
  for (cc in 1:ncol(C)) {
    xx <- freqB[,cc];
    yy <- C[,cc];
    X <- cbind(xx,yy);
    ok <- (is.finite(X) & -1 < X & X < 30);
    ok <- ok[,1] & ok[,2];
    X <- X[ok,];
    smoothScatter(X, xlim=c(0,1), ylim=Clim, xlab=xlab, ylab=Clab);
    abline(h=0:3, lty=3, lwd=1);
    abline(v=0:2/2, lty=3, lwd=1);
  
    # Plot centers
    fit <- kmeans(X, centers=centers);
    print(fit$centers);
    points(fit$centers, pch=19, col="red");
  }
  devDone();
}
