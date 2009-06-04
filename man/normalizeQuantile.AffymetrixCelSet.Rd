%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeQuantile.AffymetrixCelSet}
\alias{normalizeQuantile.AffymetrixCelSet}
\alias{AffymetrixCelSet.normalizeQuantile}
\alias{normalizeQuantile.AffymetrixCelSet}
\alias{normalizeQuantile,AffymetrixCelSet-method}


\title{Normalizes samples to have the same empirical distribution}

\description{
 Normalizes samples to have the same empirical distribution.
}

\usage{\method{normalizeQuantile}{AffymetrixCelSet}(this, path=NULL, name="normQuantile", subsetToUpdate=NULL, typesToUpdate=NULL, xTarget=NULL, subsetToAvg=subsetToUpdate, typesToAvg=typesToUpdate, ..., verbose=FALSE)}

\arguments{
  \item{path}{The path where to save the normalized data files.
    If \code{\link[base]{NULL}}, a default name is used.}
  \item{name}{The name of the normalized data set, which will also be
    part of the default path.}
  \item{subsetToUpdate}{The probes to be updated.
    If \code{\link[base]{NULL}}, all probes are updated.}
  \item{typesToUpdate}{Types of probes to be updated.}
  \item{xTarget}{A \code{\link[base]{numeric}} \code{\link[base]{vector}}.  The empirical distribution
    to which all arrays should be normalized to.}
  \item{subsetToAvg}{The probes to calculate average empirical
    distribution over.  If a single \code{\link[base]{numeric}} in (0,1), then this
    fraction of all probes will be used.
    If \code{\link[base]{NULL}}, all probes are considered.}
  \item{typesToAvg}{Types of probes to be used when calculating the
    average empirical distribution.
    If \code{"pm"} and \code{"mm"} only perfect-match and mismatch
    probes are used, respectively. If \code{"pmmm"} both types are used.
  }
  \item{...}{Additional arguments passed to \code{normalizeQuantile()}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{double}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.light]{normalizeQuantileRank.numeric}}
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
