############################################################################
# EXAMPLES:
#
# // Run all system test sets for the Test3 chip type
# Rscript testScripts/launch.R --pattern=Test3
#
# // Run a random system (default) test set
# Rscript testScripts/launch.R --order=random --nbrOfSets 1
#
# // Run a random replication test set
# Rscript testScripts/launch.R --group=replication --order=random --nbrOfSets 1
#
# // Run one of the complete analyses
# Rscript testScripts/launch.R --group=complete --pattern=GSE12702
############################################################################
library("R.utils");

# Override default settings with command line arguments  
args <- commandArgs(asValues=TRUE, excludeReserved=TRUE, excludeEnvVars=TRUE);
print(args);

printf("Hostname: %s\n", System$getHostname());

library("aroma.affymetrix");

path <- system.file(package="aroma.affymetrix");
path <- Arguments$getReadablePath(path);

path <- file.path(path, "testScripts/R");
path <- Arguments$getReadablePath(path);

pathname <- file.path(path, "launchUtils.R");
pathname <- Arguments$getReadablePathname(pathname);

source(pathname);

do.call(launchTestGroups, args);

############################################################################
# HISTORY:
# 2012-10-19
# o ROBUSTNESS: Better pathname validation.
# o Adding system details to output at the beginning.
# 2012-09-14
# o Created.
############################################################################
