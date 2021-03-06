library("aroma.affymetrix")

ovars <- ls(all.names=TRUE)
oplan <- future::plan()

## Setup dataset
dataset <- "GSE8605"
chipType <- "Mapping10K_Xba142"

csR <- AffymetrixCelSet$byName(dataset, chipType=chipType)
csR <- csR[1:2]
print(csR)

checksum <- NULL

strategies <- future:::supportedStrategies()
strategies <- setdiff(strategies, "multiprocess")
if (require("future.batchtools")) {
  strategies <- c(strategies, "batchtools_local")
  if (any(grepl("PBS_", names(Sys.getenv())))) {
    strategies <- c(strategies, "batchtools_torque")
  } else if (any(grepl("SGE_", names(Sys.getenv())))) {
    strategies <- c(strategies, "batchtools_sge")
  }
}

message("Future strategies: ", paste(sQuote(strategies), collapse = ", "))
mprint(future::sessionDetails())
mprint(list(
  availableCores = future::availableCores(which = "all"),
  availableWorkers = future::availableWorkers(which = "all")
))

for (strategy in strategies) {
  message(sprintf("*** Using %s futures ...", sQuote(strategy)))

  future::plan(strategy)
  tags <- c("*", strategy)

  ## (a) Process a single array
  bpn <- BasePositionNormalization(csR[1], target="zero", tags=c(tags, "one-array"))
  print(bpn)
  csC1 <- process(bpn, verbose=-10)
  print(csC1)
  csC1z <- getChecksumFileSet(csC1)
  print(csC1z[[1]])

  ## Compare file checksum to previous runs
  checksumT <- readChecksum(csC1z[[1]])
  if (is.null(checksum)) checksum <- checksumT
  stopifnot(identical(checksumT, checksum))


  ## (b) Process two arrays
  bpn <- BasePositionNormalization(csR, target="zero", tags=tags)
  print(bpn)
  csC <- process(bpn, verbose=-10)
  print(csC)
  csCz <- getChecksumFileSet(csC)
  print(csCz[[1]])
  res <- equals(csC1z[[1]], csCz[[1]])
  if (!res) throw(res)

  message(sprintf("*** Using %s futures ... DONE", sQuote(strategy)))
}


## CLEANUP
future::plan(oplan)
rm(list=setdiff(ls(all.names=TRUE), ovars))
