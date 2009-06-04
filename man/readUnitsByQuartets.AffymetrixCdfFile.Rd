%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.getCellQuartets.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readUnitsByQuartets.AffymetrixCdfFile}
\alias{readUnitsByQuartets.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.readUnitsByQuartets}
\alias{readUnitsByQuartets.AffymetrixCdfFile}
\alias{readUnitsByQuartets,AffymetrixCdfFile-method}


\title{Gets the cell quartets for each base position}

\description{
  Gets the cell quartets for each base position.
}

\usage{\method{readUnitsByQuartets}{AffymetrixCdfFile}(this, units=NULL, ..., verbose=FALSE)}

\arguments{
 \item{units}{Subset of units to be queried. If \code{\link[base]{NULL}}, all units are used.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{vector}} of \code{\link[base]{factor}}s.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}
\keyword{internal}
\keyword{methods}
