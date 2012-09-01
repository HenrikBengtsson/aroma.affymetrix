library("aroma.affymetrix");

# Local functions
getDataSet <- function(path, ...) {
  chipType <- basename(path);
  cdf <- tryCatch({
    AffymetrixCdfFile$byChipType(chipType);
  }, error=function(ex) {
    AffymetrixCdfFile$byChipType(chipType, tags=".*");
  });
  ds <- AffymetrixCelSet$byPath(path, cdf=cdf);
  ds;
}

getRawDataSetPath <- function(dataSet=NULL, tags=NULL, chipType=NULL, ...) {
  if (is.null(dataSet) && is.null(chipType)) return(NULL);
  dataSetF <- paste(c(dataSet, tags), collapse=",");
  path <- file.path("rawData", dataSetF, chipType);
  path <- Arguments$getWritablePath(path);
  path;
} # getRawDataSetPath()


getGeoUrl <- function(dataSet, ...) {
  sprintf("ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/series/%s/%s_RAW.tar", dataSet, dataSet);
} # getGeoUrl()


getHapMapUrlPath <- function(chipType=c("Mapping50K_Hind240", "Mapping50K_Xba240", "GenomeWideSNP_6"), ...) {
  # Argument 'chipType':
  chipType <- match.arg(chipType);

  dirs <- c("Mapping50K_Hind240"="affy100k", "Mapping50K_Xba240"="affy100k",
            "GenomeWideSNP_6"="hapmap3_affy6.0");
  dir <- dirs[chipType];

  sprintf("http://hapmap.ncbi.nlm.nih.gov/downloads/raw_data/%s", dir);
} # getHapMapUrlPath()
 

# url <- getHapMapUrl("CEU_NA06985", chipType="Mapping50K_Hind240");
getHapMapUrl <- function(sampleName, chipType=c("Mapping50K_Hind240", "Mapping50K_Xba240"), ...) {
  # Argument 'sampleName':
  sampleName <- Arguments$getCharacters(sampleName);

  # Argument 'chipType':
  chipType <- match.arg(chipType);

  tags <- c("Mapping50K_Hind240"="_HIND", "Mapping50K_Xba240"="_XBA");
  tag <- tags[chipType];

  # Create URL
  urlPath <- getHapMapUrlPath(chipType);
  sprintf("%s/%s%s.CEL.gz", urlPath, sampleName, tag);
} # getHapMapUrl()


# downloadHapMapSample(dataSet, chipType=chipType, sampleName=sampleName);
downloadHapMapSample <- function(dataSet, tags=NULL, chipType, sampleName, ..., gunzip=TRUE, skip=TRUE) {
  path <- getRawDataSetPath(dataSet=dataSet, tags=tags, chipType=chipType, ...);
  url <- getHapMapUrl(sampleName, chipType=chipType);
  filenameGZ <- basename(url);
  filename <- gsub("[.]gz$", "", filenameGZ);
  pathname <- file.path(path, filename);
  if (skip && isFile(pathname)) {
    return(pathname);
  }

  pathD <- dirname(path);
  pathnameGZ <- downloadFile(url, filename=filenameGZ, path=pathD);
  if (regexpr("[.]gz$", pathnameGZ) != -1) {
    gunzip(pathnameGZ);
  } 

  opwd <- getwd();
  on.exit(if (!is.null(opwd)) setwd(opwd));
  setwd(pathD);
  pathname <- arrangeCelFilesByChipType(filename, path=".");
  setwd(opwd); opwd <- NULL;
  pathname <- file.path(path, filename);
  pathname;
} # downloadHapMapSample()


downloadHapMapSamples <- function(..., sampleNames) {
  pathnames <- sapply(sampleNames, FUN=function(sampleName) {
    downloadHapMapSample(..., sampleName=sampleName);
  });

  path <- dirname(pathnames[1]);
  getDataSet(path);
} # downloadHapMapSamples()


downloadHapMapDataSet <- function(dataSet, tags=NULL, chipType=c("GenomeWideSNP_6"), ..., skip=TRUE) {
  # Argument 'chipType':
  chipType <- match.arg(chipType);

  population <- dataSet;
  dataSet <- paste(c("HapMap", population), collapse=",");
  path <- getRawDataSetPath(dataSet=dataSet, tags=tags, chipType=chipType, ...);
  # Already downloaded?
  ds <- getDataSet(path);
  if (skip && nbrOfFiles(ds) > 0) {
    return(ds);
  }

  urlPath <- getHapMapUrlPath(chipType);
  filenameD <- sprintf("%s.tgz", population);
  pathD <- dirname(path);
  url <- file.path(urlPath, filenameD);
  pathnameD <- downloadFile(url, filename=filenameD, path=pathD);
  untar(pathnameD, exdir=pathD);

  # Arrange by chip type
  opwd <- getwd();
  on.exit(if (!is.null(opwd)) setwd(opwd));
  setwd(pathD);
  pathnames <- list.files(pattern="[.](cel|CEL)$");
  pathnamesD <- arrangeCelFilesByChipType(pathnames, path=".");
  setwd(opwd); opwd <- NULL;

  getDataSet(path);
} # downloadHapMapDataSet()


downloadGeoRawDataSet <- function(dataSet, tags=NULL, chipType, ..., url=getGeoUrl(dataSet), gunzip=TRUE, skip=TRUE) {
  path <- getRawDataSetPath(dataSet=dataSet, tags=tags, chipType=chipType, ...);

  # Already downloaded?
  ds <- getDataSet(path);
  if (skip && nbrOfFiles(ds) > 0) {
    return(ds);
  }

  pathname <- downloadFile(url, ...);
  pathD <- dirname(path);
  untar(pathname, exdir=pathD);
  if (gunzip) {
    pathnames <- list.files(path=pathD, pattern="[.]gz$", full.names=TRUE);
    sapply(pathnames, FUN=gunzip);
  }

  # Arrange by chip type
  opwd <- getwd();
  on.exit(if (!is.null(opwd)) setwd(opwd));
  setwd(pathD);
  pathnames <- list.files(pattern="[.](cel|CEL)$");
  pathnamesD <- arrangeCelFilesByChipType(pathnames, path=".");
  setwd(opwd); opwd <- NULL;

  getDataSet(path);
} # downloadGeoRawDataSet()


############################################################################
# HISTORY:
# 2012-08-31
# o Added as a utility test script of the aroma.affymetrix package.
# 2012-08-23
# o Fixed up HapMap downloader.
# 2012-08-22
# o Added GEO and HapMap downloader.
# o Created.
############################################################################
