%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getUnitIntensities.AffymetrixCelSet}
\alias{getUnitIntensities.AffymetrixCelSet}
\alias{AffymetrixCelSet.getUnitIntensities}
\alias{getUnitIntensities.AffymetrixCelSet}
\alias{getUnitIntensities,AffymetrixCelSet-method}


\title{Gets cell signals for a subset of units and a subset of arrays}

\description{
  Gets cell signals for a subset of units and a subset of arrays.
}

\usage{\method{getUnitIntensities}{AffymetrixCelSet}(this, units=NULL, ..., force=FALSE, cache=TRUE, verbose=FALSE)}

\arguments{
 \item{units}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying units to be read.
   If \code{\link[base]{NULL}}, all units are read.}
 \item{...}{Arguments passed to the low-level function for read units,
   e.g. \code{\link[affxparser]{readCelUnits}} or \code{\link[aroma.apd]{readApdUnits}}.
   If \code{units} is not already a CDF \code{\link[base]{list}} structure, these arguments
   are also passed to \code{readUnits()} of \code{\link{AffymetrixCdfFile}}.}
 \item{force}{If \code{\link[base:logical]{TRUE}}, cached values are ignored.}
 \item{verbose}{If \code{\link[base:logical]{TRUE}}, progress details are printed, otherwise not.
   May also be a \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns a named \code{\link[base]{list}} structure.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
