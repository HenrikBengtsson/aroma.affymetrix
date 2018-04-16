###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod isSnpChip
#
# @title "Static method to check if a chip is a mapping (SNP) chip"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE if the chip type refers to a SNP array, otherwise @FALSE.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isSnpChip", "AffymetrixCdfFile", function(this, ...) {
  chipType <- getChipType(this)

  # First some hardwired return values
  if (regexpr("^Mapping(10K|50K|250K)_.*$", chipType) != -1)
    return(TRUE)

  if (regexpr("^Cent(Hind|Xba).*$", chipType) != -1)
    return(TRUE)

  if (regexpr("^GenomeWideSNP_.*$", chipType) != -1)
    return(TRUE)

  if (regexpr("^Cyto.*Array$", chipType) != -1)
    return(TRUE)

  # Then, check for genotype units
  types <- getUnitTypes(this, ...)
  hasSnpUnits <- any(types == 2)

  hasSnpUnits
}, private=TRUE)



###########################################################################/**
# @RdocMethod getSnpNames
#
# @title "Gets the names of the SNP units"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character @vector.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "getCnNames".
#   Internally, @seemethod "getUnitTypes".
#   is used.
#
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getSnpNames", "AffymetrixCdfFile", function(this, ...) {
  types <- getUnitTypes(this, ...)
  units <- (types == 2)
  getUnitNames(this, units=units, ...)
}, private=TRUE)




###########################################################################/**
# @RdocMethod getCnNames
#
# @title "Gets the names of the CN units"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character @vector.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "getSnpNames".
#   Internally, @seemethod "getUnitTypes".
#   is used.
#
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getCnNames", "AffymetrixCdfFile", function(this, ...) {
  types <- getUnitTypes(this, ...)
  units <- (types == 5)
  getUnitNames(this, units=units, ...)
}, private=TRUE)




###########################################################################/**
# @RdocMethod nbrOfSnps
#
# @title "Gets the number of SNPs"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to @seemethod "getSnpNames".}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author "HB"
#
# \seealso{
#   Internally, @seemethod "getSnpNames" is used to identify SNPs.
#
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfSnps", "AffymetrixCdfFile", function(this, ...) {
  length(getSnpNames(this, ...))
}, private=TRUE)
