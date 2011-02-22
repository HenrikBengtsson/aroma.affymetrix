# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE


.First.lib <- function(libname, pkgname) {
##   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##   # Loading/installing affxparser
##   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##   # Load 'affxparser'
##   res <- require("affxparser");
## 
##   # Not installed?
##   if (!res) {
##     if (interactive()) {
##       cat("Package 'affxparser' is not available or could not be loaded. Will now try to install it from Bioconductor (requires working internet connection):\n");
##       source("http://www.bioconductor.org/biocLite.R");
##       biocLite("affxparser");
##       # Assert that the package can be successfully loaded
##       res <- require("affxparser");
##       if (!res) {
##         throw("Package 'affxparser' could not be loaded. Please install it from Bioconductor.");
##       }
##     } else {
##       warning("Package 'affxparser' could not be loaded. Without it ", pkgname, " will not work. Please install it from Bioconductor.");
##     }
##   }


  pkg <- AromaAffymetrix(pkgname);
  assign(pkgname, pkg, pos=getPosition(pkg));

  packageStartupMessage(getName(pkg), " v", getVersion(pkg), " (", 
    getDate(pkg), ") successfully loaded. See ?", pkgname, " for help.");

  # Setup package
  .setupAromaAffymetrix(pkg);
}
