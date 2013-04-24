%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Transform.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{getOutputDataSet.Transform}
\alias{getOutputDataSet.Transform}
\alias{Transform.getOutputDataSet}
\alias{getOutputDataSet,Transform-method}

\title{Gets the transformed data set}

\description{
 Gets the transformed data set, if processed.
}

\usage{
\method{getOutputDataSet}{Transform}(this, ..., verbose=FALSE)
}

\arguments{
  \item{...}{Not used.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, any in-memory cached results are ignored.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link[aroma.core]{AromaMicroarrayDataSet}} or \code{\link[base]{NULL}}.
}

\seealso{
  For more information see \code{\link{Transform}}.
}
\keyword{internal}
\keyword{methods}
