for (zzz in 0) {

# Setup verbose output
verbose <- Arguments$getVerbose(-2)
timestampOn(verbose)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define an example dataset using this path
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Find any SNP dataset
path <- NULL
if (is.null(path))
  break

if (!exists("ds")) {
  ds <- AffymetrixCelSet$fromFiles(path)
}
print(ds)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get the genome information for this data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- getCdf(ds)
gi <- getGenomeInformation(cdf)
print(gi)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create a CN model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("model")) {
  model <- MbeiCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE)
}
print(model)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit all units on one chromosome 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chr <- "8"
units <- getUnitIndices(gi, chromosome=chr)
cat(sprintf("Fits %d SNPs on chromosome %s.\n", length(units), chr))
unitsDone <- fit(model, units=units, moreUnits=10, verbose=verbose)
nDone <- length(unitsDone)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get some information about this chromsome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chrStats <- getChromosomeStats(gi)[chr,]
print(chrStats)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get the estimated chip effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ces <-getChipEffects(model)

# Array to be investigated
ss <- grep("NA12802", getNames(ds))
ce <- getFile(ces, ss)

# Re-calculate the average chip-effect across arrays
ceAvg <- getAverageFile(ces, force=(nDone > 0), moreUnits=10, verbose=verbose)
ceAvg$mergeStrands <- ce$mergeStrands
ceAvg$combineAlleles <- ce$combineAlleles

yRef <- unlist(ceAvg[units], use.names=FALSE)
y <- unlist(ce[units], use.names=FALSE)
cn <- y/yRef	

x <- getPositions(gi, units=units, na.rm=FALSE)
# Make sure to plot in order
x <- x/10^6
xlab <- "Physical position (Mb)"
ylab <- expression(log[2](theta[i]/theta[avg]))
xlim <- c(0,chrStats["max"]/10^6);
ylim <- c(-1,1)*2
main <- sprintf("CN estimates for chromosome %s", chr)
plot(NA, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main)
abline(h=0, col="blue")
abline(h=log(1:4/2, base=2), col="blue", lty=2)
points(x, log(cn,2), pch=176)
stext(sprintf("%d SNPs", sum(!is.na(x))), side=3, pos=1, line=-1, cex=0.8)
stextChipType(cdf)

} # for (zzz in 0)
rm(zzz)
