%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenomeInformation.AFFX.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{GenomeInformation$fromCdf}
\alias{GenomeInformation$fromCdf}
\alias{fromCdf.GenomeInformation}
\alias{GenomeInformation.fromCdf}
\alias{fromCdf,GenomeInformation-method}

\title{Static method to define a genome information set from a CDF}

\description{
  Static method to define a genome information set from a CDF.
}

\usage{
## Static method (use this):
## GenomeInformation$fromCdf(cdf, ...)

## Don't use the below:
\method{fromCdf}{GenomeInformation}(static, cdf, ...)
}

\arguments{
 \item{cdf}{An \code{\link{AffymetrixCdfFile}}.}
 \item{...}{Additional argument passed to \code{\link[aroma.affymetrix:byChipType.GenomeInformation]{*byChipType}()}.}
}

\value{
  Returns a \code{\link{GenomeInformation}} object.
}

\seealso{
  \code{\link[aroma.affymetrix:byChipType.GenomeInformation]{*byChipType}()}.
  For more information see \code{\link{GenomeInformation}}.
}
\keyword{internal}
\keyword{methods}
