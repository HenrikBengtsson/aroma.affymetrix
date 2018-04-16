# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFUNCT
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEPRECATED
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TODO, but still used alot internally.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getData", "AffymetrixCelFile", function(this, ...) {
##  .Deprecated("readRawData")
  readRawData(this, ...)
}, protected=TRUE, deprecated=TRUE)


setMethodS3("nbrOfArrays", "AffymetrixCelSet", function(this, ...) {
##  .Deprecated("length")
  length(this, ...)
}, protected=TRUE)

setMethodS3("nbrOfArrays", "AffymetrixCnChpSet", function(this, ...) {
##  .Deprecated("length")
  length(this, ...)
}, protected=TRUE)

setMethodS3("nbrOfArrays", "CnagCfhSet", function(this, ...) {
##  .Deprecated("length")
  length(this, ...)
}, protected=TRUE)

setMethodS3("nbrOfArrays", "DChipDcpSet", function(this, ...) {
##  .Deprecated("length")
  length(this, ...)
}, protected=TRUE)
