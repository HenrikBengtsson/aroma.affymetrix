library("aroma.affymetrix");
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
# Spatial intensity plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ae <- ArrayExplorer(csR);
setColorMaps(ae, "sqrt,yellow");
print(ae);
stopifnot(identical(unname(getArrays(ae)), getNames(csR)));
process(ae, verbose=log);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spatial probe log-ratio plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cfR <- getAverageFile(csR, verbose=log);
reporter <- SpatialReporter(csR, reference=cfR);
ylab <- expression(log[2](y/y[R]));
for (array in 1:nbrOfArrays(csR)) {
  df <- getFile(csR, array);

  toPNG(getFullName(df), tags="spatial,rowMedians", width=300, aspectRatio=8/3, {
    plotMargins(reporter, array=array, ylim=c(-1,1)*0.2, ylab=ylab, margins="rows", rotate=90);
  });

  toPNG(getFullName(df), tags="spatial,colMedians", width=800, aspectRatio=3/8, {
    plotMargins(reporter, array=array, ylim=c(-1,1)*0.2, ylab=ylab, margins="columns", rotate=0);
  });
}

addColorMap(reporter, "log2center,rainbow");
process(reporter, zrange=c(-2,2), verbose=log);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spatial residual plots test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaPlm(csR);
print(plm);
fit(plm, verbose=log);
rs <- calculateResidualSet(plm, verbose=log);
ae <- ArrayExplorer(rs);
setColorMaps(ae, c("log2,log2neg,rainbow", "log2,log2pos,rainbow"));
print(ae);
stopifnot(identical(unname(getArrays(ae)), getNames(csR)));
process(ae, interleaved="auto", verbose=log);
