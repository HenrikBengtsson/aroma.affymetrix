%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doCRMAv1.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{doCRMAv1}
\alias{doCRMAv1.default}
\alias{doCRMAv1}
\alias{doCRMAv1.AffymetrixCelSet}
\alias{doASCRMAv1}
\alias{doASCRMAv1.default}

\title{Estimation and assessment of raw copy numbers at the single locus level (CRMA v1)}

\description{
 Estimation and assessment of raw copy numbers at the single locus level (CRMA v1) based on [1].
 The algorithm is processed in bounded memory, meaning virtually
 any number of arrays can be analyzed on also very limited computer
 systems.
}

\usage{
  \method{doCRMAv1}{AffymetrixCelSet}(csR, shift=+300, combineAlleles=TRUE, lengthRange=NULL, arrays=NULL, drop=TRUE,
  verbose=FALSE, ...)
  \method{doCRMAv1}{default}(dataSet, ..., verbose=FALSE)
  \method{doASCRMAv1}{default}(...)
}

\arguments{
 \item{csR, dataSet}{An \code{\link{AffymetrixCelSet}} (or the name of an \code{\link{AffymetrixCelSet}}).}
 \item{shift}{An tuning parameter specifying how much to shift the
  probe signals before probe summarization.}
 \item{combineAlleles}{A \code{\link[base]{logical}} specifying whether allele probe pairs
  should be summed before modeling or not.}
 \item{lengthRange}{An optional \code{\link[base]{numeric}} vector of length two passed
  to \code{\link{FragmentLengthNormalization}}.}
 \item{arrays}{A \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the subset of arrays
  to process.  If \code{\link[base]{NULL}}, all arrays are considered.}
 \item{drop}{If \code{\link[base:logical]{TRUE}}, the summaries are returned, otherwise
  a named \code{\link[base]{list}} of all intermediate and final results.}
 \item{verbose}{See \code{\link{Verbose}}.}
 \item{...}{Additional arguments used to set up \code{\link{AffymetrixCelSet}} (when argument \code{dataSet} is specified).}
}

\value{
  Returns a named \code{\link[base]{list}}, iff \code{drop == FALSE}, otherwise
  only \code{\link{ChipEffectSet}} object.
}


\section{Allele-specific or only total-SNP signals}{
  If you wish to obtain allele-specific estimates for SNPs, which
  are needed to call genotypes or infer parent-specific copy numbers,
  then use argument \code{combineAlleles=FALSE}.  Total copy number
  signals are still available.
  If you know for certain that you will not use allele-specific
  estimates, you will get slightly less noisy signals
  (very small difference) if you use \code{combineAlleles=TRUE}.

  \code{doASCRMAv1(...)} is a wrapper for
  \code{doCRMAv1(..., combineAlleles=FALSE)}.
}

\references{
 [1] H. Bengtsson, R. Irizarry, B. Carvalho & T.P. Speed.
     \emph{Estimation and assessment of raw copy numbers at the
     single locus level},
     Bioinformatics, 2008.\cr
}

\seealso{
 For CRMA v2 (recommended by authors), which is a single-array
 improvement over CRMA v1, see \code{\link{doCRMAv2}}().
}

\author{Henrik Bengtsson}

