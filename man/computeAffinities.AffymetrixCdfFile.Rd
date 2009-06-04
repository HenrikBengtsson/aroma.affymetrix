%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.computeAffinities.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{computeAffinities.AffymetrixCdfFile}
\alias{computeAffinities.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.computeAffinities}
\alias{computeAffinities.AffymetrixCdfFile}
\alias{computeAffinities,AffymetrixCdfFile-method}


\title{Calculates probe affinities from sequence}

\description{
 Calculates probe affinities from sequence.

Adapted from \code{\link[gcrma]{compute.affinities}} in the \pkg{gcrma} package.
Attempts to find the tab-separated probe sequence file associated with
a particular CDF, and matches sequence to probe index in order to assign
an affinity to each probe.
}

\usage{\method{computeAffinities}{AffymetrixCdfFile}(this, paths=NULL, safe=TRUE, force=FALSE, verbose=FALSE, ...)}

\arguments{
  \item{paths}{A \code{\link[base]{character}} variable containing location(s) to look for
   the probe sequence file.  Multiple locations should be separated by
   semi-colons.}
  \item{force}{If \code{\link[base:logical]{FALSE}}, cached results is returned, if available.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of (log2) probe affinities, of length equal
 to the total number of features on the array.
}

\author{
  Ken Simpson (ksimpson[at]wehi.edu.au).
}
\keyword{internal}
\keyword{methods}
