%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getUnitIndices.GenomeInformation}
\alias{getUnitIndices.GenomeInformation}
\alias{GenomeInformation.getUnitIndices}
\alias{getUnitIndices.GenomeInformation}
\alias{getUnitIndices,GenomeInformation-method}


\title{Gets unit indices ordered along the chromosome}

\description{
  Gets unit indices ordered along the chromosome.
}

\usage{\method{getUnitIndices}{GenomeInformation}(this, ..., na.rm=TRUE)}

\arguments{
 \item{...}{Arguments passed to \code{\link[aroma.affymetrix:getData.GenomeInformation]{*getData}()}.}
 \item{na.rm}{If \code{\link[base:logical]{TRUE}}, non-defined unit indices are excluded, otherwise
     not.}
}

\value{
  Returns an \code{\link[base]{integer}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getData.GenomeInformation]{*getData}()}.
  For more information see \code{\link{GenomeInformation}}.
}
\keyword{internal}
\keyword{methods}
