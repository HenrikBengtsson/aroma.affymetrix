%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getIntensities.AffymetrixCelSet}
\alias{getIntensities.AffymetrixCelSet}
\alias{AffymetrixCelSet.getIntensities}
\alias{getIntensities,AffymetrixCelSet-method}

\title{Gets cell intensities from a set of cells and a set of arrays}

\description{
  Gets cell intensities from a set of cells and a set of arrays.
}

\usage{
\method{getIntensities}{AffymetrixCelSet}(this, ...)
}

\arguments{
 \item{...}{Passed to \code{\link[aroma.affymetrix:getData.AffymetrixCelSet]{*getData}()}.}
}

\value{
  Returns a \code{\link[base]{numeric}} \eqn{NxK} matrix, where \eqn{N} is the number of
  cells read, and \eqn{K} is the number of arrays in the data set.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
