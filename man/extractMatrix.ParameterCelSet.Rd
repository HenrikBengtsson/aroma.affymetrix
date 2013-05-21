%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ParameterCelSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{extractMatrix.ParameterCelSet}
\alias{extractMatrix.ParameterCelSet}
\alias{ParameterCelSet.extractMatrix}
\alias{extractMatrix,ParameterCelSet-method}

\title{Extract data as a matrix for a set of arrays}

\description{
 Extract data as a matrix for a set of arrays.
}

\usage{
\method{extractMatrix}{ParameterCelSet}(this, units=NULL, ..., field=c("intensities", "stdvs", "pixels"),
  returnUgcMap=FALSE, drop=FALSE, verbose=FALSE)
}

\arguments{
  \item{units}{(The subset of units to be matched.
    If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Passed to \code{\link[base]{subset}}() operating on the UGC map.}
  \item{field}{The field to be extracted.}
  \item{returnUgcMap}{If \code{\link[base:logical]{TRUE}}, the (unit, group, cell) map is returned
    as an attribute.}
  \item{drop}{If \code{\link[base:logical]{TRUE}}, singleton dimensions are dropped.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an JxK \code{\link[base]{double}} \code{\link[base]{matrix}} where J is the number of units,
 and K is the number of arrays.
 The names of the columns are the names of the arrays.
 No names are set for the rows.
 The rows are ordered according to \code{units} arguments.
}

\author{Henrik Bengtsson}

\seealso{
  \code{\link[aroma.affymetrix:extractDataFrame.ParameterCelSet]{*extractDataFrame}()}.
  For more information see \code{\link{ParameterCelSet}}.
}
\keyword{internal}
\keyword{methods}
