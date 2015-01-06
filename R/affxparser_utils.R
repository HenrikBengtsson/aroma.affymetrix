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

