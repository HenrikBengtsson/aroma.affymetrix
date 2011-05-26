if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chipType <- "HuGene-1_0-st-v1";

footer <- list(
  createdOn = format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
  createdBy = list(
    fullname = "Henrik Bengtsson", 
    email = "hb@stat.berkeley.edu"
  ),
  srcFiles = list()
);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- AffymetrixCdfFile$byChipType(chipType, tags=".*");
print(cdf);

ptb <- AffymetrixProbeTabFile$byChipType(chipType);
print(ptb);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allocate aroma cell sequence file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acs <- AromaCellSequenceFile$allocateFromCdf(cdf, tags="*,HB20090927");
print(acs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Import data from the Affymetrix probe-tab file (contains only PMs)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
importFrom(acs, ptb, verbose=log);

# Infer MM from PM sequences?  Will give an error if no MMs
inferMmFromPm(acs, cdf=cdf, verbose=log); # DOES NOT WORK
