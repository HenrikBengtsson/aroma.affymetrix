library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");

chipType <- "MOUSEDIVm520650";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadCDF(ar, chipType);
verbose && cat(verbose, "CDF: ", pathname);

pathname <- downloadACS(ar, chipType);
verbose && cat(verbose, "ACS: ", pathname);

pathname <- downloadUGP(ar, chipType);
verbose && cat(verbose, "UGP: ", pathname);

pathname <- downloadUFL(ar, chipType);
verbose && cat(verbose, "UFL: ", pathname);

verbose && exit(verbose);
