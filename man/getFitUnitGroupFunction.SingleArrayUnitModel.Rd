%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SingleArrayUnitModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getFitUnitGroupFunction.SingleArrayUnitModel}
\alias{getFitUnitGroupFunction.SingleArrayUnitModel}
\alias{SingleArrayUnitModel.getFitUnitGroupFunction}
\alias{getFitUnitGroupFunction.SingleArrayUnitModel}
\alias{getFitUnitGroupFunction,SingleArrayUnitModel-method}


\title{Static method to get the low-level function that fits the PLM}

\description{
 Static method to get the low-level function that fits the PLM.
 Any subclass model must provide this method, which should return
 a \code{\link[base]{function}} that accepts a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length K, where K
 is the number of probes.
}

\synopsis{getFitUnitGroupFunction.SingleArrayUnitModel(...)}
\usage{SingleArrayUnitModel$getFitUnitGroupFunction()}

\arguments{
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{function}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SingleArrayUnitModel}}.
}
\keyword{internal}
\keyword{methods}
