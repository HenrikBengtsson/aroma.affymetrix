%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getPositions.GenomeInformation}
\alias{getPositions.GenomeInformation}
\alias{GenomeInformation.getPositions}
\alias{getPositions.GenomeInformation}
\alias{getPositions,GenomeInformation-method}


\title{Gets the physical positions for a set of units}

\description{
  Gets the physical positions for a set of units.
}

\usage{\method{getPositions}{GenomeInformation}(this, ..., na.rm=FALSE)}

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
  \code{\link[aroma.affymetrix:getUnitIndices.GenomeInformation]{*getUnitIndices}()}.
  For more information see \code{\link{GenomeInformation}}.
}
\keyword{internal}
\keyword{methods}
