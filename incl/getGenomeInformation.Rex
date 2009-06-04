
for (zzz in 0) {

# Define a CDF for a SNP chip
if (!exists("cdf")) {
  cdf <- AffymetrixCdfFile$fromChipType(".*Sty.*")
}
print(cdf)

# Get the genome information
gi <- getGenomeInformation(cdf)
print(gi)

# Get the units on chromosome X order position by default
# (First call will be slow be data is being cached)
units <- getUnitIndices(gi, chromosome="X")
nunits <- length(units)
cat(sprintf("Number of units on chromosome X: %d\n", nunits))
str(units)

# Plot the SNP density for this chromosome
plotDensity(gi, chromosome="X");

} # for (zzz in 0)
rm(zzz)
