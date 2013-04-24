%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{isPm.AffymetrixCdfFile}
\alias{isPm.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.isPm}
\alias{isPm,AffymetrixCdfFile-method}

\title{Checks which cells (probes) are PMs and not}

\description{
 Checks which cells (probes) are PMs and not.
}

\usage{
\method{isPm}{AffymetrixCdfFile}(this, units=NULL, force=FALSE, cache=TRUE, ..., verbose=FALSE)
}

\arguments{
  \item{units}{The units to be read. If \code{\link[base]{NULL}}, all units are read.}
  \item{...}{Additional arguments passed to \code{\link[affxparser]{readCdfUnits}}.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, cached results are ignored.}
  \item{cache}{If \code{\link[base:logical]{TRUE}}, results are cached.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{logical}} \code{\link[base]{vector}} of length K, where K equals the total number
 of cells in the requested units.  Note that the cells are ordered as
 they occur in the units, that is, \emph{not} in incremental order.
}

\section{Caching}{
  This method caches a \code{\link[base]{logical}} \code{\link[base]{vector}} of length N, when N equals the
  number of cells on the array. The size of this vector is approximately
  4*N bytes.  The vector indicates if a cell is a PM or not.
}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
