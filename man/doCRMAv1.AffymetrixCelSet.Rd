%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doCRMAv1.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{doCRMAv1.AffymetrixCelSet}
\alias{doCRMAv1.AffymetrixCelSet}
\alias{AffymetrixCelSet.doCRMAv1}
\alias{doCRMAv1,AffymetrixCelSet-method}
\alias{doCRMAv1.character}
\alias{doCRMAv1}
\alias{doASCRMAv1}

\title{Estimation and assessment of raw copy numbers at the single locus level (CRMA v1)}

\description{
 Estimation and assessment of raw copy numbers at the single locus level (CRMA v1) based on [1].
 The algorithm is processed in bounded memory, meaning virtually
 any number of arrays can be analyzed on also very limited computer
 systems.

 \code{doASCRMAv1(...)} is a wrapper for
 \code{doCRMAv1(..., combineAlleles=FALSE)}.
}

\usage{
\method{doCRMAv1}{AffymetrixCelSet}(csR, shift=+300, combineAlleles=TRUE, lengthRange=NULL, arrays=NULL, ..., drop=TRUE, ram=NULL, verbose=FALSE)
}

\arguments{
 \item{csR}{An \code{\link{AffymetrixCelSet}}.}
 \item{shift}{An tuning parameter specifying how much to shift the
  probe signals before probe summarization.}
 \item{combineAlleles}{A \code{\link[base]{logical}} specifying whether allele probe pairs
  should be summed before modelling or not.}
 \item{lengthRange}{}
 \item{arrays}{A \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the subset of arrays
  to run RMA on.  If \code{\link[base]{NULL}}, all arrays are considered.}
 \item{plm}{}
 \item{drop}{If \code{\link[base:logical]{TRUE}}, the RMA summaries are returned, otherwise
  a named \code{\link[base]{list}} of all intermediate and final results.}
 \item{verbose}{See \code{\link{Verbose}}.}
 \item{...}{Not used.}
}

\value{
  Returns a named \code{\link[base]{list}}, iff \code{drop == FALSE}, otherwise
  only \code{\link{ChipEffectSet}} object (containing the RMA summaries).
}

\references{
 [1] H. Bengtsson, R. Irizarry, B. Carvalho & T.P. Speed.
     \emph{Estimation and assessment of raw copy numbers at the
     single locus level},
     Bioinformatics, 2008.\cr
}

\seealso{
 \code{\link{doCRMAv2}}().
}

\author{Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
