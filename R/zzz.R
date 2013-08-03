# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE;

# WORKAROUND: In order for the package to work with the most recent
# version of R devel, which automatically add namespaces to packages
# who do not have one, we explicitly have specify the following.
# /HB 2011-07-27
# R.utils:
##cat <- R.utils::cat;
##getOption <- R.utils::getOption;
##lapply <- R.utils::lapply;

# R.filesets:
append <- R.filesets::append;
##sapply <- R.filesets::sapply;

# aroma.core:
apply <- aroma.core::apply;
colMeans <- aroma.core::colMeans;
colSums <- aroma.core::colSums;
library <- aroma.core::library;
require <- aroma.core::require;
.Machine <- aroma.core::.Machine;


.onAttach <- function(libname, pkgname) {
  pkg <- AromaAffymetrix(pkgname);
  assign(pkgname, pkg, pos=getPosition(pkg));

  # Setup package
  .setupAromaAffymetrix(pkg);

  startupMessage(pkg);
}
