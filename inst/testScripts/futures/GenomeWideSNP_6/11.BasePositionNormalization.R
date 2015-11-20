library("aroma.affymetrix")

ovars <- ls(all.names=TRUE)

## Setup dataset
csC <- AffymetrixCelSet$byName("GSE13372,testset", chipType="GenomeWideSNP_6,Full")
csC <- csC[1:2]
print(csC)

strategies <- c("eager", "lazy")
if (future::supportsMulticore()) strategies <- c(strategies, "multicore")

for (strategy in strategies) {
  message(sprintf("*** Using %s futures ...", sQuote(strategy)))

  future::plan(strategy)
  tags <- c("*", strategy)

  csC1 <- csC[1]
  bpn <- BasePositionNormalization(csC1, target="zero", tags=c(tags, "one-array"))
  print(bpn)
  csB1 <- process(bpn, verbose=-10)
  print(csB1)

  bpn <- BasePositionNormalization(csC, target="zero", tags=tags)
  print(bpn)
  csB <- process(bpn, verbose=-10)
  print(csB)

  ## Assert same results by comparing file checksum
  csB1z <- getChecksumFileSet(csB1)
  csBz <- getChecksumFileSet(csB)
  res <- equals(csBz[[1]], csB1z[[1]])
  if (!res) throw(res)

  message(sprintf("*** Using %s futures ... DONE", sQuote(strategy)))
}


## CLEANUP
rm(list=setdiff(ls(all.names=TRUE), ovars))
future::plan("eager")
