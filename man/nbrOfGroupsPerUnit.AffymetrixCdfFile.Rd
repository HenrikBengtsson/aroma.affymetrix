%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.COUNTS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{nbrOfGroupsPerUnit.AffymetrixCdfFile}
\alias{nbrOfGroupsPerUnit.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.nbrOfGroupsPerUnit}
\alias{nbrOfGroupsPerUnit,AffymetrixCdfFile-method}
\alias{AffymetrixCdfFile.getUnitSizes}
\alias{getUnitSizes.AffymetrixCdfFile}
\alias{getUnitSizes,AffymetrixCdfFile-method}


\title{Gets the number of groups in each unit}

\description{
 Gets the number of groups in each unit.
}

\usage{
\method{nbrOfGroupsPerUnit}{AffymetrixCdfFile}(this, units=NULL, ..., force=FALSE, verbose=FALSE)
}

\arguments{
  \item{units}{The units of interest. If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Not used.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{vector}} of \code{\link[base]{integer}}s.
}

\details{
  Once read from file, this information is cached in memory for efficiency.
  The cache can be cleared by calling \code{gc(cdf)}.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
