%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpProbeAffinityFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCellIndices.SnpProbeAffinityFile}
\alias{getCellIndices.SnpProbeAffinityFile}
\alias{SnpProbeAffinityFile.getCellIndices}
\alias{getCellIndices,SnpProbeAffinityFile-method}

\title{Retrieves tree list of cell indices for a set of units}

\description{
  Retrieves tree list of cell indices for a set of units from the associated CDF.
}

\usage{
\method{getCellIndices}{SnpProbeAffinityFile}(this, ..., unlist=FALSE)
}

\arguments{
 \item{...}{Additional arguments passed to \code{getCellIndices()}
            of \code{\link{ProbeAffinityFile}}.}
 \item{unlist}{If \code{\link[base:logical]{TRUE}}, the cell indices are returned as a \code{\link[base]{vector}}.}
}

\value{
  Returns a \code{\link[base]{list}} structure, where each element corresponds to a unit.
  If argument \code{unlist=TRUE} is passed, an \code{\link[base]{integer}} \code{\link[base]{vector}} is returned.
}

\seealso{
  For more information see \code{\link{SnpProbeAffinityFile}}.
}


\keyword{internal}
\keyword{methods}
