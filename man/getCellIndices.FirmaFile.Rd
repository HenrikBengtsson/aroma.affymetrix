%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FirmaFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCellIndices.FirmaFile}
\alias{getCellIndices.FirmaFile}
\alias{FirmaFile.getCellIndices}
\alias{getCellIndices.FirmaFile}
\alias{getCellIndices,FirmaFile-method}


\title{Retrieves tree list of cell indices for a set of units}

\description{
  Retrieves tree list of cell indices for a set of units from the associated CDF.
}

\usage{\method{getCellIndices}{FirmaFile}(this, ..., .cache=TRUE)}

\arguments{
 \item{...}{Additional arguments passed to \code{getCellIndices()}
            of \code{\link{AffymetrixCdfFile}}.}
 \item{.cache}{Ignored.}
}

\value{
  Returns a \code{\link[base]{list}} structure, where each element corresponds to a unit.
  If argument \code{unlist=TRUE} is passed, an \code{\link[base]{integer}} \code{\link[base]{vector}} is returned.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{FirmaFile}}.
}


\keyword{internal}
\keyword{methods}
