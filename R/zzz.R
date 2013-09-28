# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE;

# WORKAROUND: In order for the package to work with the most recent
# version of R devel, which automatically add namespaces to packages
# who do not have one, we explicitly have specify the following.
# /HB 2011-07-27
append <- R.filesets::append;
apply <- aroma.core::apply;
colMeans <- aroma.core::colMeans;
colSums <- aroma.core::colSums;
library <- aroma.core::library;
require <- aroma.core::require;
.Machine <- aroma.core::.Machine;


.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname);
  pkg <- AromaAffymetrix(pkgname);
  assign(pkgname, pkg, envir=ns);
} # .onLoad()


.onAttach <- function(libname, pkgname) {
  pkg <- get(pkgname, envir=getNamespace(pkgname));

  # Setup package
  .setupAromaAffymetrix(pkg);

  startupMessage(pkg);
}
