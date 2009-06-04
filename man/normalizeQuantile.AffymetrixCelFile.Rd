%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.normalizeQuantile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeQuantile.AffymetrixCelFile}
\alias{normalizeQuantile.AffymetrixCelFile}
\alias{AffymetrixCelFile.normalizeQuantile}
\alias{normalizeQuantile.AffymetrixCelFile}
\alias{normalizeQuantile,AffymetrixCelFile-method}


\title{Normalizes the probe intensities to a target empirical distribution}

\description{
 Normalizes the probe intensities to a target empirical distribution.
}

\usage{\method{normalizeQuantile}{AffymetrixCelFile}(this, path=file.path("normQuantile", getChipType(this)), xTarget, subsetToUpdate=NULL, typesToUpdate=NULL, ..., overwrite=FALSE, skip=!overwrite, verbose=FALSE)}

\arguments{
  \item{path}{The path where to save the normalized data files.}
  \item{xTarget}{A \code{\link[base]{numeric}} \code{\link[base]{vector}}.  The empirical distribution
    to which all arrays should be normalized to.}
  \item{subsetToUpdate}{The indices of the probes to be updated.
    If \code{\link[base]{NULL}}, all are updated.}
  \item{typesToUpdate}{Types of probes to be updated.  For more details,
    see argument \code{types} of \code{identifyCells()} for the
    \code{\link{AffymetrixCdfFile}} class.}
  \item{...}{Additional arguments passed to \code{normalizeQuantile()}.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already normalized arrays are overwritten,
    unless skipped, otherwise an error is thrown.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the array is not normalized if it already exists.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the normalized \code{\link{AffymetrixCelFile}} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.light]{normalizeQuantile.numeric}}
  For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
