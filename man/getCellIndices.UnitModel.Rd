%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  UnitModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCellIndices.UnitModel}
\alias{getCellIndices.UnitModel}
\alias{UnitModel.getCellIndices}
\alias{getCellIndices,UnitModel-method}

\title{Gets the cell indices unit by unit}

\description{
 Gets the cell indices unit by unit for all or a subset of units (probesets).
}

\usage{
\method{getCellIndices}{UnitModel}(this, ..., verbose=FALSE)
}

\arguments{
  \item{...}{Additional arguments passed to \code{getCellIndices()}
    of the \code{\link{AffymetrixCdfFile}} of the input data set.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the \code{\link[base]{list}} structure consisting of CDF cell indices.
}

\seealso{
  For more information see \code{\link{UnitModel}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
