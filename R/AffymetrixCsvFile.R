setConstructorS3("AffymetrixCsvFile", function(..., sep=",", .verify=TRUE) {
  this <- extend(AffymetrixTabularFile(..., .verify=FALSE), "AffymetrixCsvFile");

  if (.verify)
    verify(this, ...);
  this;
})


setMethodS3("getExtensionPattern", "AffymetrixCsvFile", function(static, ...) {
  "[.](csv|CSV)$";
}, static=TRUE, protected=TRUE)


setMethodS3("findByChipType", "AffymetrixCsvFile", function(static, chipType, pattern=sprintf("^%s.*[.](csv|CSV)$", chipType), ...) {
  findByChipType.AffymetrixTabularFile(static, chipType=chipType, pattern=pattern, ...);
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
