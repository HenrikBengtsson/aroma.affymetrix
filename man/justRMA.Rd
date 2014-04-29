%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  justRMA.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{justRMA}
\alias{justRMA}

\alias{justRMA.AffymetrixCelSet}
\alias{justRMA.default}

\title{Robust Multichip Analysis (RMA) reproducing the affy package}

\description{
 Robust Multichip Analysis (RMA) reproducing the affy package as far as possible.  The RMA method is described in [1].

 The algorithm is processed in bounded memory, meaning a very large
 number of arrays can be analyzed on also very limited computer systems.
 The only limitation is the amount of memory required to load the final
 chip-effect estimates into memory (as a \code{\link[Biobase]{ExpressionSet}}).
}

\usage{
  \method{justRMA}{AffymetrixCelSet}(csR, flavor=c("oligo", "affyPLM"), ..., verbose=FALSE)
  \method{justRMA}{default}(...)
}

\arguments{
 \item{csR}{An \code{\link{AffymetrixCelSet}}.}
 \item{flavor}{A \code{\link[base]{character}} string specifying the estimators used in
  the RMA summarization step.}
 \item{...}{Additional arguments passed to \code{\link{doRMA}}() used internally.}
 \item{verbose}{See \code{\link{Verbose}}.}
}

\value{
  Returns an annotated \code{\link[Biobase]{ExpressionSet}}.
}

\section{Reproducibility of affy}{
  This implementation of the RMA method reproduces \code{\link[affy]{justRMA}}
  in \pkg{affy} package quite well.  It does so by still using a
  constant memory profile, i.e. it is possible to use this implementation
  to run RMA on a much large data set than what is possible with
  \pkg{affy}.  At least 20-50 \emph{times} more samples should be
  doable, if not more.
}

\references{
 [1] Irizarry et al.
     \emph{Summaries of Affymetrix GeneChip probe level data}.
     NAR, 2003, 31, e15.\cr
}

\seealso{
  \code{\link{doRMA}}().
}

\author{Henrik Bengtsson}
\keyword{methods}
