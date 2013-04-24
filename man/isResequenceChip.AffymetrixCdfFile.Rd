%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.getCellQuartets.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{isResequenceChip.AffymetrixCdfFile}
\alias{isResequenceChip.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.isResequenceChip}
\alias{isResequenceChip,AffymetrixCdfFile-method}

\title{Static method to check if a CDF is for a resequencing (CustomSeq) chip}

\description{
  Static method to check if a CDF is for a resequencing (CustomSeq) chip.  Note, this method is not bullet proof.  Several
  resequencing CDF does not carry that information.  For such, we add
  tests based on their chip type, as we become aware of them.
}

\usage{
\method{isResequenceChip}{AffymetrixCdfFile}(this, ...)
}

\arguments{
 \item{...}{Not used.}
}

\value{
  Returns \code{\link[base:logical]{TRUE}} if the chip type refers to a resequence array,
  otherwise \code{\link[base:logical]{FALSE}}.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
