%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doCRMAv2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{doCRMAv2.AffymetrixCelSet}
\alias{doCRMAv2.AffymetrixCelSet}
\alias{AffymetrixCelSet.doCRMAv2}
\alias{doCRMAv2,AffymetrixCelSet-method}
\alias{doCRMAv2.character}
\alias{doCRMAv2}
\alias{doASCRMAv2}

\title{A single-array preprocessing method for estimating full-resolution raw copy numbers from all Affymetrix genotyping arrays (CRMA v2)}

\description{
 A single-array preprocessing method for estimating full-resolution raw copy numbers from all Affymetrix genotyping arrays (CRMA v2) based on [1].
 The algorithm is processed in bounded memory, meaning virtually
 any number of arrays can be analyzed on also very limited computer
 systems.

 \code{doASCRMAv2(...)} is a wrapper for
 \code{doCRMAv2(..., combineAlleles=FALSE)}.
}

\usage{
\method{doCRMAv2}{AffymetrixCelSet}(csR, combineAlleles=TRUE, lengthRange=NULL, arrays=NULL, ..., plm=c("AvgCnPlm", "RmaCnPlm"), drop=TRUE, ram=NULL, verbose=FALSE)
}

\arguments{
 \item{csR}{An \code{\link{AffymetrixCelSet}}.}
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
 [1] H. Bengtsson, P. Wirapati & T.P. Speed.
     \emph{A single-array preprocessing method for estimating
     full-resolution raw copy numbers from all Affymetrix
     genotyping arrays including GenomeWideSNP 5 \& 6},
     Bioinformatics, 2009.\cr
}

\seealso{
 \code{\link{doCRMAv1}}().
}

\author{Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
