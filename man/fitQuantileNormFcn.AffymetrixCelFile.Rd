%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.fitQuantileNormFcn.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitQuantileNormFcn.AffymetrixCelFile}
\alias{fitQuantileNormFcn.AffymetrixCelFile}
\alias{AffymetrixCelFile.fitQuantileNormFcn}
\alias{fitQuantileNormFcn,AffymetrixCelFile-method}

\title{Fits quantile normalization functions for the arrays in the data set}

\description{
 Fits quantile normalization functions for the arrays in the data set.
}

\usage{
\method{fitQuantileNormFcn}{AffymetrixCelFile}(this, yTarget, subset=NULL, ..., spar=NULL, nknots=1024,
  verbose=FALSE)
}

\arguments{
  \item{yTarget}{The target probe signals.}
  \item{subset}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} specifying the indices of the
     subset of probes to be used to fit the normalization function.}
  \item{spar, nknots}{Control parameters passed to
     \code{\link[stats]{smooth.spline}}.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a normalization \code{\link[base]{function}}.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
