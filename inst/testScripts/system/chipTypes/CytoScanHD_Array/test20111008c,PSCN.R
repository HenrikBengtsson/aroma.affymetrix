library("aroma.affymetrix");

verbose <- Arguments$getVerbose(-50, timestamp=TRUE);

dataSetName <- "Affymetrix_2011-CytoScanHD";
tags <- "ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY";
chipType <- "CytoScanHD_Array";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- list();
dsList$total <- AromaUnitTotalCnBinarySet$byName(dataSetName, tags=tags, chipType=chipType);
dsList$fracB <- AromaUnitFracBCnBinarySet$byName(dataSetName, tags=tags, chipType=chipType);
print(dsList);

dfR <- getAverageFile(dsList$total);
thetaR <- extractMatrix(dfR, drop=TRUE);

array <- 7;
dfList <- lapply(dsList, FUN=getFile, array);
sampleName <- getFullName(dfList$total);

data <- sapply(dfList, FUN=extractMatrix);
data[,"total"] <- 2 * data[,"total"] / thetaR;
str(data);

ugp <- getAromaUgpFile(dsList$total);

chr <- 8;
units <- getUnitsOnChromosome(ugp, chr);
x <- ugp[units,2,drop=TRUE];
x <- x/1e6;

dataT <- data[units,];

C <- dataT[,"total"];
B <- dataT[,"fracB"];

chrTag <- sprintf("Chr%02d", chr);
toPNG(sampleName, tags=c(chrTag, "TCN+BAF"), width=1024, aspectRatio=0.5, {
  par(mar=c(2.2,3,1.2,1));
  subplots(2, ncol=1);
  plot(x, C, pch=".", ylim=c(0,6));
  stext(side=3, pos=0, sampleName);
  stext(side=3, pos=1, chrTag);
  plot(x, B, pch=".", ylim=c(0,1));
})
