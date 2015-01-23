library("aroma.affymetrix")

if (setupExampleData(aroma.affymetrix, dataset="FusionSDK_Test3", mustWork=FALSE)) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- AffymetrixCdfFile$byChipType("Test3")
  print(cdf)

  # Unit names, types, ...
  names <- getUnitNames(cdf)
  str(names)

  units <- indexOf(cdf, names=names[42:40])
  stopifnot(all.equal(units, 42:40))

  types <- getUnitTypes(cdf)
  print(table(types))

  data <- readUnits(cdf, units=1:10)
  str(data)

  data <- readDataFrame(cdf, units=1:10)
  str(data)

  md5 <- getChecksumFile(cdf)
  print(md5)

  cdfM <- getMonocellCdf(cdf, verbose=-100)
  print(cdfM)
} # if (... "AffymetrixDataTestFiles")
