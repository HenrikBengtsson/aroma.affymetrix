%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.getAlleleProbePairs.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{getAlleleProbePairs.AffymetrixCdfFile}
\alias{getAlleleProbePairs.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.getAlleleProbePairs}
\alias{getAlleleProbePairs,AffymetrixCdfFile-method}

\title{Gets the indices of probepairs with the same pair of SNP nucleotides}

\description{
  Gets the indices of probepairs with the same pair of SNP nucleotides.
  Note that the order of allele A and allele B is irrelevant.
  For instance, all probepairs with nucleotides (A,T) are calibrated
  together with all probepairs with nucleotides (T,A) reversed.
}

\usage{
\method{getAlleleProbePairs}{AffymetrixCdfFile}(this, units=NULL, ignoreOrder=TRUE, force=FALSE, verbose=FALSE, ...)
}

\arguments{
 \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
 \item{...}{Not used.}
}

\value{
  Returns a named \code{\link[base]{list}} where each element is a two-column \code{\link[base]{matrix}} where
  the column names are the nucleotides for the two alleles.
}

\section{Benchmarking}{
  On an IBM Thinkpad A31 with 1.8GHz and 1GB RAM:
  \itemize{
   \item{Mapping10K\_Xba142}{10208 units & 432964 cells: 11 seconds.}
   \item{Mapping50K\_Xba240}{58960 SNPs & 589,600 (PMA,PMB) probe pairs: 11 seconds.}
  }
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}
\keyword{internal}
\keyword{methods}
