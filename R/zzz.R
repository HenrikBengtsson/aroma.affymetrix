# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE


.First.lib <- function(libname, pkgname) {
  pkg <- AromaAffymetrix(pkgname);
  assign(pkgname, pkg, pos=getPosition(pkg));

  packageStartupMessage(getName(pkg), " v", getVersion(pkg), " (", 
    getDate(pkg), ") successfully loaded. See ?", pkgname, " for help.");

  # Setup package
  .setupAromaAffymetrix(pkg);
}
