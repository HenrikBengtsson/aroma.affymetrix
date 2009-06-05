%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readCdfDataFrame.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{readCdfDataFrame}
\alias{readCdfDataFrame}


 \title{Reads units (probesets) from an Affymetrix CDF file}

 \usage{readCdfDataFrame(filename, units=NULL, groups=NULL, cells=NULL, fields=NULL, drop=TRUE, verbose=0)}

 \description{
  Reads units (probesets) from an Affymetrix CDF file. Gets all or a subset of units (probesets).
 }

 \arguments{
  \item{filename}{The filename of the CDF file.}
  \item{units}{An \code{\link[base]{integer}} \code{\link[base]{vector}} of unit indices
    specifying which units to be read.  If \code{\link[base]{NULL}}, all are read.}
  \item{groups}{An \code{\link[base]{integer}} \code{\link[base]{vector}} of group indices
    specifying which groups to be read.  If \code{\link[base]{NULL}}, all are read.}
  \item{cells}{An \code{\link[base]{integer}} \code{\link[base]{vector}} of cell indices
    specifying which cells to be read.  If \code{\link[base]{NULL}}, all are read.}
  \item{fields}{A \code{\link[base]{character}} \code{\link[base]{vector}} specifying what fields to read.
    If \code{\link[base]{NULL}}, all unit, group and cell fields are returned.}
  \item{drop}{If \code{\link[base:logical]{TRUE}} and only one field is read, then a \code{\link[base]{vector}}
    (rather than a single-column \code{\link[base]{data.frame}}) is returned.}
  \item{verbose}{An \code{\link[base]{integer}} specifying the verbose level. If 0, the
    file is parsed quietly.  The higher numbers, the more details.}
 }

 \value{
   An NxK \code{\link[base]{data.frame}} or a \code{\link[base]{vector}} of length N.
 }

 \author{Henrik Bengtsson}

 \examples{
##############################################################
if (require("AffymetrixDataTestFiles")) {            # START #
##############################################################

# Find any CDF file
cdfFile <- NULL
tryCatch({
  cdfFile <- findCdf()
}, error = function(ex) {})

if (!is.null(cdfFile)) {
  units <- 101:120
  fields <- c("unit", "unitName", "group", "groupName", "cell")
  df <- readCdfDataFrame(cdfFile, units=units, fields=fields)
  stopifnot(identical(sort(unique(df$unit)), units))
 
  fields <- c("unit", "unitName", "unitType")
  fields <- c(fields, "group", "groupName")
  fields <- c(fields, "x", "y", "cell", "pbase", "tbase")
  df <- readCdfDataFrame(cdfFile, units=units, fields=fields)
  stopifnot(identical(sort(unique(df$unit)), units))
}

##############################################################
}                                                     # STOP #
##############################################################

}

 \seealso{
   For retrieving the CDF as a \code{\link[base]{list}} structure, see
   \code{\link[affxparser]{readCdfUnits}}.
 }

 \references{
   [1] Affymetrix Inc, Affymetrix GCOS 1.x compatible file formats,
       June 14, 2005.
       \url{http://www.affymetrix.com/support/developer/}
 }




\keyword{file}
\keyword{IO}
\keyword{internal}
