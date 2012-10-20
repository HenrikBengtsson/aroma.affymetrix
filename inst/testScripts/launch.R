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
path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "launchUtils.R");
source(pathname);

# Override default settings with command line arguments  
args <- commandArgs(asValues=TRUE, excludeReserved=TRUE, excludeEnvVars=TRUE);
print(args);

do.call(launchTestGroups, args);

############################################################################
# HISTORY:
# 2012-09-14
# o Created.
############################################################################
