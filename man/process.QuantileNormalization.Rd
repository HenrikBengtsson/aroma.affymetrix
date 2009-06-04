%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QuantileNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{process.QuantileNormalization}
\alias{process.QuantileNormalization}
\alias{QuantileNormalization.process}
\alias{process.QuantileNormalization}
\alias{process,QuantileNormalization-method}


\title{Normalizes the data set}

\description{
 Normalizes the data set.
}

\usage{\method{process}{QuantileNormalization}(this, ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{...}{Not used.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, data already normalized is re-normalized,
      otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{double}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{QuantileNormalization}}.
}
\keyword{internal}
\keyword{methods}
