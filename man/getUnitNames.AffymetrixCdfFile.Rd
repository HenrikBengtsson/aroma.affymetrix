%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getUnitNames.AffymetrixCdfFile}
\alias{getUnitNames.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.getUnitNames}
\alias{getUnitNames.AffymetrixCdfFile}
\alias{getUnitNames,AffymetrixCdfFile-method}


\title{Gets the names of each unit}

\description{
 Gets the names of each unit.
}

\usage{\method{getUnitNames}{AffymetrixCdfFile}(this, units=NULL, ...)}

\arguments{
  \item{units}{The units of interest. If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{vector}} of \code{\link[base]{character}} names.
}

\details{
  Once read from file, this information is cached in memory for efficiency.
  The cache can be cleared by calling \code{gc(cdf)}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
