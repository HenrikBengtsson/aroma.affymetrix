%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCellIndices.ChipEffectFile}
\alias{getCellIndices.ChipEffectFile}
\alias{ChipEffectFile.getCellIndices}
\alias{getCellIndices,ChipEffectFile-method}

\title{Retrieves tree list of cell indices for a set of units}

\description{
  Retrieves tree list of cell indices for a set of units from the associated CDF.
}

\usage{
\method{getCellIndices}{ChipEffectFile}(this, ..., .cache=TRUE)
}

\arguments{
 \item{...}{Additional arguments passed to \code{getCellIndices()}
            of \code{\link{AffymetrixCdfFile}}.}
 \item{.cache}{Ignored.}
}

\value{
  Returns a \code{\link[base]{list}} structure, where each element corresponds to a unit.
  If argument \code{unlist=TRUE} is passed, an \code{\link[base]{integer}} \code{\link[base]{vector}} is returned.
}

\seealso{
  For more information see \code{\link{ChipEffectFile}}.
}


\keyword{internal}
\keyword{methods}
