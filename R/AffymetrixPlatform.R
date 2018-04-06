# @author "HB"
setConstructorS3("AffymetrixPlatform", function(...) {
  extend(AromaPlatform(), "AffymetrixPlatform")
})

setMethodS3("getName", "AffymetrixPlatform", function(static, ...) {
  "Affymetrix"
})

setMethodS3("getUnitNamesFile", "AffymetrixPlatform", function(static, ...) {
  AffymetrixCdfFile$byName(...)
}, static=TRUE)

setMethodS3("getUnitTypesFile", "AffymetrixPlatform", function(static, ...) {
  AffymetrixCdfFile$byName(...)
}, static=TRUE)
