##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
subset <- c(2,5,3);

ds <- doCRMAv2(dataSet, chipType=chipType, arrays=subset, verbose=log);
print(ds);

columns <- c("unitName", "chromosome", "position", "*");
df <- getFile(ds, 2);
dfTxt <- writeDataFrame(df, columns=columns, overwrite=TRUE, verbose=log);
print(dfTxt);
data <- readDataFrame(dfTxt, row=1:100);
str(data);

dfTxt <- writeDataFrame(ds, overwrite=TRUE, verbose=log);
print(dfTxt);
data <- readDataFrame(dfTxt, row=1:100);
str(data);

columns <- c("*");
dfTxt <- writeDataFrame(ds, columns=columns, overwrite=TRUE, verbose=log);
print(dfTxt);
data <- readDataFrame(dfTxt, row=1:100);
str(data);

columns <- c("unitName", "chromosome", "position", "*");
dfTxt <- writeDataFrame(ds, columns=columns, overwrite=TRUE, verbose=log);
print(dfTxt);
data <- readDataFrame(dfTxt, row=1:100);
str(data);
