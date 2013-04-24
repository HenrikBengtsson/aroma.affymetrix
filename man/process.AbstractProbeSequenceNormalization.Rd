%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractProbeSequenceNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{process.AbstractProbeSequenceNormalization}
\alias{process.AbstractProbeSequenceNormalization}
\alias{AbstractProbeSequenceNormalization.process}
\alias{process,AbstractProbeSequenceNormalization-method}

\title{Normalizes the data set}

\description{
 Normalizes the data set.
}

\usage{
\method{process}{AbstractProbeSequenceNormalization}(this, ..., ram=NULL, force=FALSE, verbose=FALSE)
}

\arguments{
  \item{...}{Not used.}
  \item{ram}{A positive \code{\link[base]{double}} scale factor specifying how much more
    memory to use relative to the default.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, data already normalized is re-normalized,
      otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{double}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson, James Bullard, Kasper Hansen, Pierre Neuvial, Elizabeth Purdom, Mark Robinson, Ken Simpson}

\seealso{
  For more information see \code{\link{AbstractProbeSequenceNormalization}}.
}
\keyword{internal}
\keyword{methods}
