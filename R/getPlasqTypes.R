getPlasqTypes <- function(cdf, ...) {
  cdf <- applyCdfGroups(cdf, cdfAddPlasqTypes);
  cdf <- applyCdfGroups(cdf, cdfGetFields, c("plasqType"));
#  cdf <- applyCdfGroups(cdf, cdfMergeAlleles);
  cdf <- base::lapply(cdf, FUN=.subset2, "groups");
  cdf;
}

############################################################################
# HISTORY:
# 2006-12-30
# o Created.
############################################################################  
