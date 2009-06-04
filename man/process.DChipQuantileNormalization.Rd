%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipQuantileNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{process.DChipQuantileNormalization}
\alias{process.DChipQuantileNormalization}
\alias{DChipQuantileNormalization.process}
\alias{process.DChipQuantileNormalization}
\alias{process,DChipQuantileNormalization-method}


\title{Normalizes the data set}

\description{
 Normalizes the data set.
}

\usage{\method{process}{DChipQuantileNormalization}(this, ..., force=FALSE, skip=TRUE, verbose=FALSE)}

\arguments{
  \item{...}{Arguments passed to
      \code{\link[aroma.light]{normalizeQuantileSpline.numeric}}.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, data already normalized is re-normalized,
      otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{double}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{DChipQuantileNormalization}}.
}
\keyword{internal}
\keyword{methods}
