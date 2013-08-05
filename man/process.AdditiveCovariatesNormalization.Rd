%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AdditiveCovariatesNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{process.AdditiveCovariatesNormalization}
\alias{process.AdditiveCovariatesNormalization}
\alias{AdditiveCovariatesNormalization.process}
\alias{process,AdditiveCovariatesNormalization-method}

\title{Normalizes the data set}

\description{
 Normalizes the data set.
}

\usage{
\method{process}{AdditiveCovariatesNormalization}(this, ..., force=FALSE, verbose=FALSE)
}

\arguments{
  \item{...}{Additional arguments passed to
    \code{\link[aroma.light]{normalizeFragmentLength}} (only for advanced users).}
  \item{force}{If \code{\link[base:logical]{TRUE}}, data already normalized is re-normalized,
      otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{double}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AdditiveCovariatesNormalization}}.
}
\keyword{internal}
\keyword{methods}
