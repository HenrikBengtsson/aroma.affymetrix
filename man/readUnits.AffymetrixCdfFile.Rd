%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readUnits.AffymetrixCdfFile}
\alias{readUnits.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.readUnits}
\alias{readUnits,AffymetrixCdfFile-method}

\title{Reads CDF data unit by unit}

\description{
 Reads CDF data unit by unit for all or a subset of units (probeset).
}

\usage{
\method{readUnits}{AffymetrixCdfFile}(this, units=NULL, ..., verbose=FALSE)
}

\arguments{
  \item{units}{The units to be read. If \code{\link[base]{NULL}}, all units are read.}
  \item{...}{Additional arguments passed to \code{\link[affxparser]{readCdfUnits}}.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the \code{\link[base]{list}} structure that \code{\link[affxparser]{readCdfUnits}} returns
 (possibly restructured).
}

\section{Caching}{
  CDF data is neither cached in memory nor on file by this method.
}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
