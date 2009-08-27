%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getUnitTypes.AffymetrixCdfFile}
\alias{getUnitTypes.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.getUnitTypes}
\alias{getUnitTypes.AffymetrixCdfFile}
\alias{getUnitTypes,AffymetrixCdfFile-method}


\title{Gets the types of a set of units}

\description{
 Gets the types of a set of units.
}

\usage{\method{getUnitTypes}{AffymetrixCdfFile}(this, units=NULL, ..., force=FALSE, .cache=TRUE, verbose=FALSE)}

\arguments{
  \item{units}{The units of interest. If \code{\link[base]{NULL}}, all units are considered.}
  \item{map}{A \code{\link[base]{character}} string specifying the mapping used.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{vector}} of \code{\link[base]{integer}}s.
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
