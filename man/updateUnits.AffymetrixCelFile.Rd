%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{updateUnits.AffymetrixCelFile}
\alias{updateUnits.AffymetrixCelFile}
\alias{AffymetrixCelFile.updateUnits}
\alias{updateUnits,AffymetrixCelFile-method}

\title{Updates CEL data unit by unit}

\description{
 Updates CEL data unit by unit for all or a subset of units.
}

\usage{
\method{updateUnits}{AffymetrixCelFile}(this, data, ...)
}

\arguments{
  \item{data}{A \code{\link[base]{list}} structure consisting of grouped data similar to
     what \code{\link[aroma.affymetrix:readUnits.AffymetrixCelFile]{*readUnits}()} returns.}
  \item{...}{Not used.}
}

\value{
 Returns the \code{\link[base]{list}} structure that \code{\link[affxparser]{updateCelUnits}} returns.
}

\author{Henrik Bengtsson}

\seealso{
  \code{\link[aroma.affymetrix:updateUnits.AffymetrixCelFile]{*updateUnits}()}.
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
