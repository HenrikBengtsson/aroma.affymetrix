%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FragmentLengthNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{process.FragmentLengthNormalization}
\alias{process.FragmentLengthNormalization}
\alias{FragmentLengthNormalization.process}
\alias{process,FragmentLengthNormalization-method}

\title{Normalizes the data set}

\description{
 Normalizes the data set.
}

\usage{
\method{process}{FragmentLengthNormalization}(this, ..., force=FALSE, verbose=FALSE)
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

\seealso{
  For more information see \code{\link{FragmentLengthNormalization}}.
}
\keyword{internal}
\keyword{methods}
