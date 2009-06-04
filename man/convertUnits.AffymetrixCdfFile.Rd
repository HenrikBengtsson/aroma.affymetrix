%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{convertUnits.AffymetrixCdfFile}
\alias{convertUnits.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.convertUnits}
\alias{convertUnits.AffymetrixCdfFile}
\alias{convertUnits,AffymetrixCdfFile-method}


\title{Gets and validates unit indices}

\description{
 Gets and validates unit indices either by unit names or by a unit indices (validation).
}

\usage{\method{convertUnits}{AffymetrixCdfFile}(this, units=NULL, keepNULL=FALSE, ...)}

\arguments{
  \item{units}{Either a \code{\link[base]{character}} \code{\link[base]{vector}} with unit names, or an \code{\link[base]{integer}}
    \code{\link[base]{vector}} with unit indices to be validated.
    If \code{\link[base]{NULL}}, all unit indices are returned.}
  \item{keepNULL}{If \code{\link[base:logical]{TRUE}}, \code{\link[base]{NULL}} returns \code{\link[base]{NULL}}.}
  \item{...}{Not used.}
}

\value{
 Returns an \code{\link[base]{integer}} \code{\link[base]{vector}} with unit indices.
 If some units are non existing, an error is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}
\keyword{internal}
\keyword{methods}
