%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpChipEffectFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCellIndices.SnpChipEffectFile}
\alias{getCellIndices.SnpChipEffectFile}
\alias{SnpChipEffectFile.getCellIndices}
\alias{getCellIndices.SnpChipEffectFile}
\alias{getCellIndices,SnpChipEffectFile-method}


\title{Retrieves tree list of cell indices for a set of units}

\description{
  Retrieves tree list of cell indices for a set of units.
}

\usage{\method{getCellIndices}{SnpChipEffectFile}(this, units=NULL, ..., unlist=FALSE, force=FALSE, .cache=TRUE, verbose=FALSE)}

\arguments{
 \item{units}{Unit indices to be retrieved.
              If \code{\link[base]{NULL}}, all units are considered.}
 \item{...}{Additional arguments passed to \code{getCellIndices()}
            of \code{\link{ChipEffectFile}}.}
 \item{unlist}{If \code{\link[base:logical]{TRUE}}, the cell indices are returned as a \code{\link[base]{vector}}.}
 \item{force}{If \code{\link[base:logical]{TRUE}}, the cell indices are re-read regardless whether
    they are already cached in memory or not.}
 \item{.cache}{(internal) If \code{\link[base:logical]{TRUE}}, the result is cached in memory.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link[base]{list}} structure, where each element corresponds to a unit.
  If argument \code{unlist=TRUE} is passed, an \code{\link[base]{integer}} \code{\link[base]{vector}} is returned.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SnpChipEffectFile}}.
}


\keyword{internal}
\keyword{methods}
