%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.BG.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{calculateParametersGsb.AffymetrixCelSet}
\alias{calculateParametersGsb.AffymetrixCelSet}
\alias{AffymetrixCelSet.calculateParametersGsb}
\alias{calculateParametersGsb,AffymetrixCelSet-method}

\title{Computes parameters for adjustment of specific binding}

\description{
 Computes parameters for adjustment of specific binding.
}

\usage{
\method{calculateParametersGsb}{AffymetrixCelSet}(this, nbrOfPms=25000, affinities=NULL, seed=NULL, ...,
  verbose=FALSE)
}

\arguments{
  \item{nbrOfPms}{The number of random PMs to use in estimation.}
  \item{affinities}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of probe affinities.}
  \item{seed}{An (optional) \code{\link[base]{integer}} specifying a temporary random seed
    to be used during processing.  The random seed is set to its original
    state when done.  If \code{\link[base]{NULL}}, it is not set.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\details{
  This method is not constant in memory! /HB 2007-03-26
}

\author{Ken Simpson}
\keyword{internal}
\keyword{methods}
