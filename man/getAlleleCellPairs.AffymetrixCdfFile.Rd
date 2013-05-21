%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.getAlleleCellPairs.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{getAlleleCellPairs.AffymetrixCdfFile}
\alias{getAlleleCellPairs.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.getAlleleCellPairs}
\alias{getAlleleCellPairs,AffymetrixCdfFile-method}

\title{Gets the cell indices of allele pairs}

\description{
  Gets the cell indices of allele pairs
  in units of type "genotyping" with 2 or 4 unit groups.
}

\usage{
\method{getAlleleCellPairs}{AffymetrixCdfFile}(this, units=NULL, stratifyBy=c("pm", "pmmm", "mm"), force=FALSE, ...,
  verbose=FALSE)
}

\arguments{
 \item{units}{A \code{\link[base]{integer}} \code{\link[base]{vector}} of units to query.
   If \code{\link[base]{NULL}}, all units are considered.}
 \item{stratifyBy}{A \code{\link[base]{character}} string specifying what type of probes
   to return.}
 \item{...}{Not used.}
 \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns a Nx2 \code{\link[base]{integer}} \code{\link[base]{matrix}} of cell indices, where each row is
  a (PMA, PMB) probe pair.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}
\keyword{internal}
\keyword{methods}
