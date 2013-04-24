%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SnpInformation$fromCdf}
\alias{SnpInformation$fromCdf}
\alias{fromCdf.SnpInformation}
\alias{SnpInformation.fromCdf}
\alias{fromCdf,SnpInformation-method}

\title{Static method to define a genome information set from a CDF}

\description{
  Static method to define a genome information set from a CDF.
}

\usage{
## Static method (use this):
## SnpInformation$fromCdf(cdf, ...)

## Don't use the below:
\method{fromCdf}{SnpInformation}(static, cdf, ...)
}

\arguments{
 \item{cdf}{An \code{\link{AffymetrixCdfFile}}.}
 \item{...}{Additional argument passed to \code{\link[aroma.affymetrix:byChipType.SnpInformation]{*byChipType}()}.}
}

\value{
  Returns a \code{\link{SnpInformation}} object.
}

\seealso{
  \code{\link[aroma.affymetrix:byChipType.SnpInformation]{*byChipType}()}.
  For more information see \code{\link{SnpInformation}}.
}
\keyword{internal}
\keyword{methods}
