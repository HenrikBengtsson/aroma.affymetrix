setMethodS3("as.AromaMicroarrayDataSetTuple", "AromaMicroarrayDataSetTuple", function(this, ...) {
  # Nothing to do
  this;
})

setMethodS3("as.AromaMicroarrayDataSetTuple", "AromaMicroarrayDataSet", function(this, ...) {
  classNames <- class(this);
  classNames <- sprintf("%sTuple", classNames);
  methodNames <- sprintf("as.%s", classNames);

  # Try first as.Nnn() method available
  for (kk in seq(along=methodNames)) {
    methodName <- methodNames[kk];
    if (exists(methodName, mode="function")) {
      fcn <- get(methodName, mode="function");
      res <- fcn(this, ...);
      return(res);
    }
  } # for (kk ...)

  throw("Failed to coerce ", class(this)[1], ", to an AromaMicroarrayDataSetTuple: ", getFullName(this));
})


##############################################################################
# HISTORY:
# 2009-11-18
# o Created.
##############################################################################

