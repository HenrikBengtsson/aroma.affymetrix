setMethodS3("allocateFromCdf", "AromaUnitGcContentFile", function(static, cdf, ...) {
  types <- "double";
  sizes <- 4;

  # NextMethod() not supported here.
  res <- allocateFromCdf.AromaUnitTabularBinaryFile(static, cdf=cdf, types=types, sizes=sizes, ...);
  res[,1] <- as.double(NA);

  res;
}, static=TRUE)



setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUnitGcContentFile", function(this, csv, rows=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csv':
  csv <- Arguments$getInstanceOf(csv, "AffymetrixNetAffxCsvFile");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing (unit name, GC content) data from ", class(csv)[1]);

  # Query CDF
  cdf <- getCdf(this);

  # Read data
  data <- readDataFrame(csv, colClassPattern=c("^(probeSetID|%GC)$"="character")); 
  unitNames <- data[,1];
  verbose && str(verbose, unitNames);
  values <- as.double(data[,2]);
  verbose && str(verbose, values);
  rm(data);
  
  # Keep the units in the CDF
  units <- indexOf(cdf, names=unitNames);
  verbose && cat(verbose, "Unit indices:");
  verbose && str(verbose, units);
  rm(unitNames);

  keep <- which(!is.na(units));
  verbose && printf(verbose, "Keeping %d of %d (%.2f%%)\n", 
       length(keep), length(units), 100*length(keep)/length(units));
  units <- units[keep];
  rm(keep);
  if (length(units) == 0) {
    warning("None of the unit names in the CSV match the ones in the CDF ('", getPathname(cdf), "'). Is the correct file ('", getPathname(csv), "'), being imported?");
  }

  # Update
  this[units,1] <- values;
  rm(values);

  verbose && exit(verbose);

  invisible(units);
})


############################################################################
# HISTORY:
# 2009-03-22
# o Created from AromaUflFile.R.
############################################################################
