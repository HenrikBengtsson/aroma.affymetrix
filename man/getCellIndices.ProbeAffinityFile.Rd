%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeAffinityFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCellIndices.ProbeAffinityFile}
\alias{getCellIndices.ProbeAffinityFile}
\alias{ProbeAffinityFile.getCellIndices}
\alias{getCellIndices,ProbeAffinityFile-method}

\title{Retrieves tree list of cell indices for a set of units}

\description{
  Retrieves tree list of cell indices for a set of units from the associated CDF.
}

\usage{
\method{getCellIndices}{ProbeAffinityFile}(this, ..., verbose=FALSE)
}

\arguments{
 \item{...}{Additional arguments passed to \code{getCellIndices()}
            of \code{\link{AffymetrixCdfFile}}.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link[base]{list}} structure, where each element corresponds to a unit.
  If argument \code{unlist=TRUE} is passed, an \code{\link[base]{integer}} \code{\link[base]{vector}} is returned.
}

\seealso{
  For more information see \code{\link{ProbeAffinityFile}}.
}


\keyword{internal}
\keyword{methods}
