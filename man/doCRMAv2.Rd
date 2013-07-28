%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doCRMAv2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{doCRMAv2}
\alias{doCRMAv2.default}
\alias{doCRMAv2}
\alias{doCRMAv2.AffymetrixCelSet}
\alias{doASCRMAv2}
\alias{doASCRMAv2.default}

\title{A single-array preprocessing method for estimating full-resolution raw copy numbers from all Affymetrix genotyping arrays (CRMA v2)}

\description{
 A single-array preprocessing method for estimating full-resolution raw copy numbers from all Affymetrix genotyping arrays (CRMA v2) based on [1].
 The algorithm is processed in bounded memory, meaning virtually
 any number of arrays can be analyzed on also very limited computer
 systems.

 We recommend CRMA v2 for estimating allele-specific as well total
 SNP signals from Affymetrix SNP chips.
}

\usage{
  \method{doCRMAv2}{AffymetrixCelSet}(csR, combineAlleles=TRUE, lengthRange=NULL, arrays=NULL,
  plm=c("AvgCnPlm", "RmaCnPlm"), drop=TRUE, verbose=FALSE, ...)
  \method{doCRMAv2}{default}(dataSet, ..., verbose=FALSE)
  \method{doASCRMAv2}{default}(...)
}

\arguments{
 \item{csR, dataSet}{An \code{\link{AffymetrixCelSet}} (or the name of an \code{\link{AffymetrixCelSet}}).}
 \item{combineAlleles}{A \code{\link[base]{logical}} specifying whether allele probe pairs
  should be summed before modelling or not.}
 \item{lengthRange}{An optional \code{\link[base]{numeric}} vector of length two passed
  to \code{\link{FragmentLengthNormalization}}.}
 \item{arrays}{A \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the subset of arrays
  to run RMA on.  If \code{\link[base]{NULL}}, all arrays are considered.}
 \item{plm}{A \code{\link[base]{character}} string specifying which type of
  probe-summarization model to used.}
 \item{drop}{If \code{\link[base:logical]{TRUE}}, the RMA summaries are returned, otherwise
  a named \code{\link[base]{list}} of all intermediate and final results.}
 \item{verbose}{See \code{\link{Verbose}}.}
 \item{...}{Additional arguments used to set up \code{\link{AffymetrixCelSet}} (when argument \code{dataSet} is specified).}
}

\value{
  Returns a named \code{\link[base]{list}}, iff \code{drop == FALSE}, otherwise
  only \code{\link{ChipEffectSet}} object (containing the RMA summaries).
}

\section{Allele-specific or only total-SNP signals}{
  If you wish to obtain allele-specific estimates for SNPs, which
  are needed to call genotypes or infer parent-specific copy numbers,
  then use argument \code{combineAlleles=FALSE}.  Total copy number
  signals are still available.
  If you know for certain that you will not use allele-specific
  estimates, you will get slightly less noisy signals
  (very small difference) if you use \code{combineAlleles=TRUE}.

  \code{doASCRMAv2(...)} is a wrapper for
  \code{doCRMAv2(..., combineAlleles=FALSE)}.
}

\references{
 [1] H. Bengtsson, P. Wirapati & T.P. Speed.
     \emph{A single-array preprocessing method for estimating
     full-resolution raw copy numbers from all Affymetrix
     genotyping arrays including GenomeWideSNP 5 \& 6},
     Bioinformatics, 2009.\cr
}

\seealso{
 For CRMA v1, which is a multi-array methods that precedes CRMA v2,
 see \code{\link{doCRMAv1}}().
}

\author{Henrik Bengtsson}

