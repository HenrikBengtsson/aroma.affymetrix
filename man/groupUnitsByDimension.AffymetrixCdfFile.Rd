%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.groupUnitsByDimension.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{groupUnitsByDimension.AffymetrixCdfFile}
\alias{groupUnitsByDimension.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.groupUnitsByDimension}
\alias{groupUnitsByDimension,AffymetrixCdfFile-method}

\title{Groups units by dimensions}

\description{
 Groups units by dimensions.
}

\usage{
\method{groupUnitsByDimension}{AffymetrixCdfFile}(this, units=NULL, ..., sort=TRUE, verbose=FALSE)
}

\arguments{
  \item{units}{An optional \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the units to be
    queried.}
  \item{...}{Not used.}
  \item{sort}{If \code{\link[base:logical]{TRUE}}, the sets are ordered by number of groups per
    units and then by the number of cells per group.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a named \code{\link[base]{list}}.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}
\keyword{internal}
\keyword{methods}
