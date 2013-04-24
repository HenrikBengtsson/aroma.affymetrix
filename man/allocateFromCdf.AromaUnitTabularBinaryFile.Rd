%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaUnitTabularBinaryFile.AFFX.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{AromaUnitTabularBinaryFile$allocateFromCdf}
\alias{AromaUnitTabularBinaryFile$allocateFromCdf}
\alias{allocateFromCdf.AromaUnitTabularBinaryFile}
\alias{AromaUnitTabularBinaryFile.allocateFromCdf}
\alias{allocateFromCdf,AromaUnitTabularBinaryFile-method}

\title{Creates an AromaUnitTabularBinaryFile mapping to a given CDF}

\description{
 Creates an AromaUnitTabularBinaryFile mapping to a given CDF.
}

\usage{
## Static method (use this):
## AromaUnitTabularBinaryFile$allocateFromCdf(cdf, ...)

## Don't use the below:
\method{allocateFromCdf}{AromaUnitTabularBinaryFile}(static, cdf, ...)
}

\arguments{
  \item{cdf}{The \code{\link{AffymetrixCdfFile}} used as a template and from
     which the (full) chip type is taken.}
  \item{path}{The path where to store the new file.}
  \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of optional tags appended to
     the filename.}
  \item{footer}{A nested named \code{\link[base]{list}} structure of additional attributes
     that are saved in the file footer after the mandatory ones.}
  \item{...}{Additional arguments passed to \code{allocate()} of
     \code{\link[aroma.core]{AromaTabularBinaryFile}}.}
}

\value{
 Returns a \code{\link[aroma.core]{AromaUnitTabularBinaryFile}} object.
}

\author{Henrik Bengtsson}

\seealso{
  To update to file footer afterwards, see \code{writeFooter()}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
