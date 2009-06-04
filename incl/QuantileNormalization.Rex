for (zzz in 0) {

# Setup verbose output
verbose <- Arguments$getVerbose(-2)
timestampOn(verbose)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define an example dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Find any dataset
path <- NULL
if (is.null(path))
  break

ds <- AffymetrixCelSet$fromFiles(path)
print(ds)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
norm <- QuantileNormalization(ds, subsetToAvg=1/3)
dsQN <- process(norm, verbose=verbose)
print(dsQN)

} # for (zzz in 0)
rm(zzz)
