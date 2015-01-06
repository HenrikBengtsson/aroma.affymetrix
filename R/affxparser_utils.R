# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Imports from affxparser
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.cdfAddPlasqTypes <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::cdfAddPlasqTypes(...)
}

.readBpmap <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readBpmap(...)
}

.readBpmapHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readBpmapHeader(...)
}

.readCdfGroupNames <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfGroupNames(...)
}

.readCdfHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfHeader(...)
}

.readCdfIsPm <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfIsPm(...)
}

.readCdf <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdf(...)
}

.readCdfQc <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfQc(...)
}

.readCdfCellIndices <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfCellIndices(...)
}

.readCdfUnits <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfUnits(...)
}

.readCdfUnitNames <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfUnitNames(...)
}

.readCdfNbrOfCellsPerUnitGroup <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfNbrOfCellsPerUnitGroup(...)
}

.createCel <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::createCel(...)
}

.compareCdfs <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::compareCdfs(...)
}

.cdfHeaderToCelHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::cdfHeaderToCelHeader(...)
}

.findCdf <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::findCdf(...)
}

.convertCdf <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::convertCdf(...)
}

.convertCel <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::convertCel(...)
}

.updateCel <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::updateCel(...)
}

.updateCelUnits <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::updateCelUnits(...)
}

.readCelRectangle <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCelRectangle(...)
}

.readPgfEnv <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readPgfEnv(...)
}

.readPgfHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readPgfHeader(...)
}

.readCcgHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCcgHeader(...)
}

.readCcg <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCcg(...)
}

.writeCdfHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::writeCdfHeader(...)
}

.writeCdfUnits <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::writeCdfUnits
(...)
}

.writeCdfQcUnits <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::writeCdfQcUnits(...)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Imports from aroma.light
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.normalizeFragmentLength <- function(...) {
  requireNamespace("aroma.light") || throw("Package not loaded: aroma.light")
  aroma.light::normalizeFragmentLength(...)
}

.normalizeQuantile <- function(...) {
  requireNamespace("aroma.light") || throw("Package not loaded: aroma.light")
  aroma.light::normalizeQuantile(...)
}

.normalizeQuantileSpline <- function(...) {
  requireNamespace("aroma.light") || throw("Package not loaded: aroma.light")
  aroma.light::normalizeQuantileSpline(...)
}

.calibrateMultiscan <- function(...) {
  requireNamespace("aroma.light") || throw("Package not loaded: aroma.light")
  aroma.light::calibrateMultiscan(...)
}

.robustSmoothSpline <- function(...) {
  requireNamespace("aroma.light") || throw("Package not loaded: aroma.light")
  aroma.light::robustSmoothSpline(...)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Imports for
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



