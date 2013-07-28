%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doRMA.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{doRMA}
\alias{doRMA.default}
\alias{doRMA}
\alias{doRMA.AffymetrixCelSet}

\title{Robust Multichip Analysis (RMA)}

\description{
 Robust Multichip Analysis (RMA) based on [1].
 The algorithm is processed in bounded memory, meaning virtually
 any number of arrays can be analyzed on also very limited computer
 systems.
 The method replicates the results of \code{\link[affyPLM]{fitPLM}}
 (package \pkg{affyPLM}) with great precision.
}

\usage{
  \method{doRMA}{AffymetrixCelSet}(csR, arrays=NULL, flavor=c("affyPLM", "oligo"), uniquePlm=FALSE, drop=TRUE,
  verbose=FALSE, ...)
  \method{doRMA}{default}(dataSet, ..., verbose=FALSE)
}

\arguments{
 \item{csR, dataSet}{An \code{\link{AffymetrixCelSet}} (or the name of an \code{\link{AffymetrixCelSet}}).}
 \item{arrays}{A \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the subset of arrays
  to run RMA on.  If \code{\link[base]{NULL}}, all arrays are considered.}
 \item{flavor}{A character string specifying what model fitting algorithm to be used, cf. \code{\link{RmaPlm}}.}
 \item{uniquePlm}{If \code{\link[base:logical]{TRUE}}, the log-additive probe-summarization model
  is done on probeset with \emph{unique} sets of probes.
  If \code{\link[base:logical]{FALSE}}, the summarization is done on "as-is" probesets as
  specified by the CDF.}
 \item{drop}{If \code{\link[base:logical]{TRUE}}, the RMA summaries are returned, otherwise
  a named \code{\link[base]{list}} of all intermediate and final results.}
 \item{verbose}{See \code{\link{Verbose}}.}
 \item{...}{Additional arguments used to set up \code{\link{AffymetrixCelSet}} (when argument \code{dataSet} is specified).}
}

\value{
  Returns a named \code{\link[base]{list}}, iff \code{drop == FALSE}, otherwise
  only \code{\link{ChipEffectSet}} object (containing the RMA summaries).
}

\references{
 [1] Irizarry et al.
     \emph{Summaries of Affymetrix GeneChip probe level data}.
     NAR, 2003, 31, e15.\cr
}

\author{Henrik Bengtsson}

