###########################################################################/**
# @RdocClass ResidualSet
#
# @title "The ResidualSet class"
#
# \description{
#  @classhierarchy
#
#  This class represents probe-level residuals from probe-level models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixCelSet".}
#   \item{probeModel}{The specific type of model, e.g. \code{"pm"}.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#   An object of this class is typically obtained through the
#   \code{getResidualSet()} method for the @see "ProbeLevelModel" class.
# }
#
#*/###########################################################################
setConstructorS3("ResidualSet", function(..., probeModel=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probeModel':
  probeModel <- match.arg(probeModel);

  extend(AffymetrixCelSet(...), "ResidualSet",
    "cached:.firstCells" = NULL,
    probeModel = probeModel
  )
})

setMethodS3("clearCache", "ResidualSet", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".firstCells")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)


setMethodS3("as.character", "ResidualSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod(generic="as.character", object=this, ...);
  params <- paste(getParametersAsString(this), collapse=", ");
  s <- c(s, sprintf("Parameters: (%s)", params));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getParameters", "ResidualSet", function(this, ...) {
  rf <- getFile(this, 1);
  getParameters(rf, ...);
})

setMethodS3("getParametersAsString", "ResidualSet", function(this, ...) {
  params <- getParameters(this);
  params <- trim(capture.output(str(params)))[-1];
  params <- gsub("^[$][ ]*", "", params);
  params <- gsub(" [ ]*", " ", params);
  params <- gsub("[ ]*:", ":", params);
  params;
}, private=TRUE)


setMethodS3("getResidualFileClass", "ResidualSet", function(static, ...) {
  ResidualFile;
}, static=TRUE, private=TRUE)


setMethodS3("byPath", "ResidualSet", function(static, ..., pattern=",residuals[.](c|C)(e|E)(l|L)$", cdf=NULL, fileClass=NULL) {
  # Argument 'cdf':
  if (!is.null(cdf)) {
    if (!inherits(cdf, "AffymetrixCdfFile"))
      throw("Argument 'cdf' must be an AffymetrixCdfFile: ", class(cdf)[1]);     
  }

  # Argument 'fileClass':
  if (is.null(fileClass))
    fileClass <- gsub("Set$", "File", class(static)[1]);

  res <- byPath.AffymetrixFileSet(static, ..., pattern=pattern, fileClass=fileClass);

  # Set CDF?
  if (!is.null(cdf))
    setCdf(res, cdf);

  res;
}, protected=TRUE, static=TRUE)


setMethodS3("fromDataSet", "ResidualSet", function(static, dataSet, path, fullname=getFullName(dataSet), cdf=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get the ResidualFile class specific for this set
  clazz <- getResidualFileClass(static);
  
  verbose && cat(verbose, "ResidualFile class: ", getName(clazz));

  verbose && enter(verbose, "Retrieving probe-level residuals from data set");
  rs <- vector("list", length(dataSet));
  verbose && cat(verbose, "Data set: ", fullname);
  for (kk in seq(dataSet)) {
    df <- getFile(dataSet, kk);
    verbose && enter(verbose, 
                           sprintf("Retrieving residual file #%d of %d (%s)",
                                               kk, length(rs), getName(df)));
    verbose && cat(verbose, "Data file class: ", class(df)[1]);
    rf <- clazz$fromDataFile(df, path=path, name=fullname, cdf=cdf, ..., 
                                                       verbose=less(verbose));
    verbose && cat(verbose, "Residual file class: ", class(rf)[1]);
    # Assert correctness
    stopifnot(inherits(rf, "ResidualFile"));

    if (is.null(cdf)) {
      verbose && enter(verbose, "Retrieving the CDF for the residual file");
      cdf <- getCdf(rf);
      verbose && exit(verbose);
    }
    rs[[kk]] <- rf;
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Create an ResidualSet
  newInstance(static, rs);
})


setMethodS3("getCellIndices", "ResidualSet", function(this, ...) {
  # Use the first residual file to get the CDF structure.
  # Note: Ideally we want to define a special CDF class doing this
  # instead of letting the data file do this. /HB 2006-12-18
  rf <- getFile(this, 1);
  getCellIndices(rf, ...);
})


setMethodS3("readUnits", "ResidualSet", function(this, units=NULL, cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading residuals unit by unit for ", nbrOfArrays(this), " arrays");

  if (is.null(cdf)) {
    verbose && enter(verbose, "Getting cell indices from CDF");
    cdf <- getCellIndices(this, units=units, ..., verbose=less(verbose));
    verbose && exit(verbose);
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  verbose && enter(verbose, "Calling readUnits() in superclass");
  res <- NextMethod("readUnits", this, units=cdf, ..., verbose=less(verbose));
  verbose && exit(verbose);

  # Get first residual file and use that to decode the read structure
  # This takes some time for a large number of units /HB 2006-10-04
  rf <- getFile(this, 1);
  res <- decode(rf, res, verbose=less(verbose));

  verbose && exit(verbose);

  res;
})


setMethodS3("updateUnits", "ResidualSet", function(this, units=NULL, cdf=NULL, data, ..., verbose=FALSE) {
  # Argument 'verbose': 
  verbose <- Arguments$getVerbose(verbose);

  # Get the CDF structure for all residual files
  if (is.null(cdf))
    cdf <- getCellIndices(this, units=units);

  # Update each file one by one
  n <- length(this);
  verbose && enter(verbose, "Updating ", n, " residual files");
  names <- getNames(this);
  verbose <- less(verbose);
  for (kk in seq(this)) {
    verbose && enter(verbose, sprintf("Array #%d of %d: %s", kk, n, names[kk]));
    rf <- getFile(this, kk);

    verbose <- less(verbose, 50);
    verbose && enter(verbose, "Extracting estimates");  # 3-4s
    dataOne <- base::lapply(data, FUN=base::lapply, function(group) {
      # eps = group$eps[,kk] = ...
      list(
        eps=.subset(.subset2(group, "eps"), kk)
      );
    });
    verbose && exit(verbose);

    verbose && enter(verbose, "Updating file");  # 6-7s ~98% in encode()
    updateUnits(rf, cdf=cdf, data=dataOne, verbose=less(verbose, 50));
    verbose && exit(verbose);
    verbose <- more(verbose, 50);

    verbose && exit(verbose);
  } # for (kk ...)
  verbose <- more(verbose);
  verbose && exit(verbose);
}, protected=TRUE)



setMethodS3("getAverageFile", "ResidualSet", function(this, ..., verbose=FALSE, indices="remaining") {
  # Argument 'indices':
  if (identical(indices, "remaining")) {
  } else if (is.null(indices)) {
    # Update only cells which stores values
    indices <- getCellIndices(this, verbose=verbose);
    indices <- unlist(indices, use.names=FALSE);
  }

  NextMethod(generic="getAverageFile", object=this, ..., indices=indices, verbose=verbose);
})


setMethodS3("findUnitsTodo", "ResidualSet", function(this, ...) {
  # Look into the last residual file since that is updated last
  rf <- getFile(this, length(this));
  findUnitsTodo(rf, ...);
})

############################################################################
# HISTORY:
# 2008-05-08
# o Made fromFiles() protected.
# 2007-12-08
# o Added argument 'cdf' to fromFiles() of ResidualSet.
# 2007-02-12
# o Created from ChipEffectSet.R.
############################################################################
