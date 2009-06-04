%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SmoothMultiarrayModel.fit.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fit.SmoothMultiarrayModel}
\alias{fit.SmoothMultiarrayModel}
\alias{SmoothMultiarrayModel.fit}
\alias{fit.SmoothMultiarrayModel}
\alias{fit,SmoothMultiarrayModel-method}


\title{Fits the model for one chromosome across samples}

\description{
 Fits the model for one chromosome across samples.
}

\usage{\method{fit}{SmoothMultiarrayModel}(this, chromosome, force=FALSE, ..., verbose=FALSE)}

\arguments{
  \item{data}{A \code{\link[base]{data.frame}} with columns \code{M} (log-ratio) and
     \code{x} (locus position).
  }
  \item{chromosome}{An \code{\link[base]{integer}} specifying the index of the chromosome to
     be fitted.}
  \item{...}{Additional arguments passed down to the internal fit function.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SmoothMultiarrayModel}}.
}
\keyword{internal}
\keyword{methods}
