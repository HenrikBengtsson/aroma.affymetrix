%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.extractMatrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{extractMatrix.AffymetrixCelSet}
\alias{extractMatrix.AffymetrixCelSet}
\alias{AffymetrixCelSet.extractMatrix}
\alias{extractMatrix,AffymetrixCelSet-method}

\title{Extract data as a matrix for a set of arrays}

\description{
 Extract data as a matrix for a set of arrays.
}

\usage{
\method{extractMatrix}{AffymetrixCelSet}(this, cells=NULL, ..., field=c("intensities", "stdvs", "pixels"),
  drop=FALSE, verbose=FALSE)
}

\arguments{
  \item{cells}{(The subset of cells to be matched.
    If \code{\link[base]{NULL}}, all cells are considered.}
  \item{...}{Not used.}
  \item{field}{The field to be extracted.}
  \item{drop}{If \code{\link[base:logical]{TRUE}}, singleton dimensions are dropped.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an JxK \code{\link[base]{double}} \code{\link[base]{matrix}} where J is the number of units,
 and K is the number of arrays.
 The names of the columns are the names of the arrays.
 No names are set for the rows.
 The rows are ordered according to \code{cells} argument.
}

\author{Henrik Bengtsson, Mark Robinson}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
