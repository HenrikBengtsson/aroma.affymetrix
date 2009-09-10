%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SingleArrayUnitModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readUnits.SingleArrayUnitModel}
\alias{readUnits.SingleArrayUnitModel}
\alias{SingleArrayUnitModel.readUnits}
\alias{readUnits.SingleArrayUnitModel}
\alias{readUnits,SingleArrayUnitModel-method}


\title{Reads data unit by unit}

\description{
 Reads data unit by unit for all or a subset of units (probeset)
 specially structured for this PLM.
}

\usage{\method{readUnits}{SingleArrayUnitModel}(this, array, units=NULL, ..., verbose=FALSE)}

\arguments{
  \item{units}{The units to be read. If \code{\link[base]{NULL}}, all units are read.}
  \item{...}{Arguments passed to \code{getCellIndices()} of the
    \code{\link{AffymetrixCdfFile}} class (if \code{cdf} was not specified),
    but also to the \code{readUnits()} method of the
    \code{\link{AffymetrixCelFile}} class.}
}

\value{
 Returns the \code{\link[base]{list}} structure that \code{readUnits()} of
 \code{\link{AffymetrixCelFile}} returns.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SingleArrayUnitModel}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
