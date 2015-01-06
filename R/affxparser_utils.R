.convertCdf <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::convertCdf(...)
} # .convertCdf()

.convertCel <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::convertCel(...)
} # .convertCel()

.updateCel <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::updateCel(...)
} # .updateCel()

.updateCelUnits <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::updateCelUnits(...)
} # .updateCelUnits()
