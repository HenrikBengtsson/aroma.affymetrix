%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AffymetrixCdfFile$findByChipType}
\alias{AffymetrixCdfFile$findByChipType}
\alias{findByChipType.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.findByChipType}
\alias{findByChipType,AffymetrixCdfFile-method}

\title{Locates a CDF file from its chip type}

\description{
 Locates a CDF file from its chip type.
}

\usage{
## Static method (use this):
## AffymetrixCdfFile$findByChipType(chipType, tags=NULL, pattern=NULL, ..., .useAffxparser=TRUE)

## Don't use the below:
\method{findByChipType}{AffymetrixCdfFile}(static, chipType, tags=NULL, pattern=NULL, ..., .useAffxparser=TRUE)
}

\arguments{
 \item{chipType}{A \code{\link[base]{character}} string.}
 \item{tags}{An optional \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
 \item{pattern}{An optional \code{\link[base]{character}} string.}
 \item{...}{Not used.}
 \item{.useAffxparser}{If \code{\link[base:logical]{TRUE}}, \code{\link[affxparser]{findCdf}} is used if
   the CDF could not be located.}
}

\value{
 Returns a pathname as a \code{\link[base]{character}} string to the first CDF file found.
 If non CDF with requested chip type was found, \code{\link[base]{NULL}} is returned.
}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
