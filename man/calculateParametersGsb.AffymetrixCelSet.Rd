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
\alias{calculateParametersGsb.AffymetrixCelSet}
\alias{calculateParametersGsb,AffymetrixCelSet-method}


\title{Computes parameters for adjustment of specific binding}

\description{
 Computes parameters for adjustment of specific binding.
}

\usage{\method{calculateParametersGsb}{AffymetrixCelSet}(this, nbrOfPms=25000, affinities=NULL, path=NULL, ..., verbose=FALSE)}

\arguments{
  \item{nbrOfPms}{The number of random PMs to use in estimation.}
  \item{affinities}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of probe affinities.}
  \item{path}{If an affinities vector is not specified,
     gives the path to a file storing the affinities.}
}

\details{
  This method is not constant in memory! /HB 2007-03-26
}

\author{
  Ken Simpson (ksimpson[at]wehi.edu.au).
}
\keyword{internal}
\keyword{methods}
