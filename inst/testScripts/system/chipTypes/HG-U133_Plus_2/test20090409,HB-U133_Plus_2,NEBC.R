library("aroma.affymetrix");

log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSet <- "Affymetrix-HeartBrain";
chipType <- "HG-U133_Plus_2";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);

bg <- RmaBackgroundCorrection(csR);
print(bg);
csB1 <- process(bg, verbose=log);
print(csB1);


# Alternative, which gives identical results
bg <- NormExpBackgroundCorrection(csR);
print(bg);
csB2 <- process(bg, verbose=log);
print(csB2);

Y1 <- extractMatrix(csB1);
Y2 <- extractMatrix(csB2);
stopifnot(identical(Y1, Y2));


# Same model, but fitted with a maximum-likelihood estimator
bg <- NormExpBackgroundCorrection(csR, method="mle");
print(bg);
csB3 <- process(bg, verbose=log);
print(csB3);

Y1 <- extractMatrix(csB1);
Y3 <- extractMatrix(csB3);
#stopifnot(identical(Y1, Y3));


