%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipSnpInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{byChipType.DChipSnpInformation}
\alias{byChipType.DChipSnpInformation}
\alias{DChipSnpInformation.byChipType}
\alias{byChipType,DChipSnpInformation-method}

\title{Defines a DChipSnpInformation object by chip type}

\description{
 Defines a DChipSnpInformation object by chip type.
}

\usage{
\method{byChipType}{DChipSnpInformation}(static, chipType, version=NULL, ...)
}

\arguments{
 \item{chipType}{A \code{\link[base]{character}} string.}
 \item{version}{An optional \code{\link[base]{character}} string specifying the version
   string, if more than one version is available.}
 \item{...}{Additional arguments passed to \code{\link[aroma.affymetrix:findByChipType.DChipSnpInformation]{*findByChipType}()}.}
}

\value{
 Returns an \code{\link{DChipSnpInformation}} object.
 If no file was not found, an error is thrown.
}

\seealso{
  For more information see \code{\link{DChipSnpInformation}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
