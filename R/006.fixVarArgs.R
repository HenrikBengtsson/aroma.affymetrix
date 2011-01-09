# Added '...' to some base functions. These will later be
# turned into default functions by setMethodS3().

write <- appendVarArgs(write);
getPackageName <- appendVarArgs(getPackageName);

# So that one can call plotRle(qa, show.names=FALSE)
boxplot.stats <- appendVarArgs(boxplot.stats);

writeCdf <- appendVarArgs(writeCdf);

############################################################################
# HISTORY:
# 2011-01-09 [HB]
# o Added appendVarArgs(writeCdf).
# 2010-02-08 [HB]
# o Added appendVarArgs(boxplot.stats) so that one can pass argument
#   'show.names' to bxp() via plotRle().
# 2007-02-27 [HB]
# o BUG FIX: Removed explicit reference to 'base' etc again. The reason is 
#   that if a previous package already modified, say, write(), to become a 
#   generic function, that was overwritten again when this package was 
#   loaded.
# 2007-02-23 [KS]
# o Make explicit reference to 'base' - this is safer, in case of colMeans()
#   defined by user or other packages.
# 2006-03-24
# o Created to please R CMD check.
############################################################################
