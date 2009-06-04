for (zzz in 0) {

# Find any dataset
path <- NULL
if (is.null(path))
  break

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define a dataset object based on all CEL files in a directory
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- AffymetrixCelSet$fromFiles(path)
print(ds)

# Keep at most three arrays for this example
ds <- extract(ds, 1:min(3,nbrOfArrays(ds)))
print(ds)

} # for (zzz in 0)
rm(zzz)
