%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doGCRMA.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{doGCRMA}
\alias{doGCRMA.default}
\alias{doGCRMA}
\alias{doGCRMA.AffymetrixCelSet}

\title{Robust Multichip Analysis (GCRMA)}

\description{
 Robust Multichip Analysis (GCRMA) based on [1].
 The algorithm is processed in bounded memory, meaning virtually
 any number of arrays can be analyzed on also very limited computer
 systems.
 The method replicates the results of \code{\link[gcrma]{gcrma}}
 (package \pkg{gcrma}) with great precision.
}

\usage{
  \method{doGCRMA}{AffymetrixCelSet}(csR, arrays=NULL, type=c("fullmodel", "affinities"), uniquePlm=FALSE, drop=TRUE, verbose=FALSE, ...)
  \method{doGCRMA}{default}(dataSet, ..., verbose=FALSE)
}

\arguments{
 \item{csR, dataSet}{An \code{\link{AffymetrixCelSet}} (or the name of an \code{\link{AffymetrixCelSet}}).}
 \item{arrays}{A \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the subset of arrays
  to run RMA on.  If \code{\link[base]{NULL}}, all arrays are considered.}
 \item{type}{A \code{\link[base]{character}} string specifying what type of model to
  use for the GCRMA background correction.
  For more details, see \code{\link{GcRmaBackgroundCorrection}}.}
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
 [1] Z. Wu, R. Irizarry, R. Gentleman, F.M. Murillo & F. Spencer.
     \emph{A Model Based Background Adjustment for Oligonucleotide
     Expression Arrays}, JASA, 2004.\cr
}

\author{Henrik Bengtsson}

