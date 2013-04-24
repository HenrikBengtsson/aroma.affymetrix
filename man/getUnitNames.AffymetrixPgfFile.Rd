%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixPgfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getUnitNames.AffymetrixPgfFile}
\alias{getUnitNames.AffymetrixPgfFile}
\alias{AffymetrixPgfFile.getUnitNames}
\alias{getUnitNames,AffymetrixPgfFile-method}

\title{Gets the names of each unit}

\description{
 Gets the names of each unit.
}

\usage{
\method{getUnitNames}{AffymetrixPgfFile}(this, units=NULL, ...)
}

\arguments{
  \item{units}{The units of interest. If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{vector}} of \code{\link[base]{character}} names.
}

\details{
  Once read from file, this information is cached in memory for efficiency.
  The cache can be cleared by calling \code{gc(pgf)}.
}

\seealso{
  For more information see \code{\link{AffymetrixPgfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
