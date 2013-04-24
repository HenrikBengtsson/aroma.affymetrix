%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultiArrayUnitModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{MultiArrayUnitModel$getFitUnitGroupFunction}
\alias{MultiArrayUnitModel$getFitUnitGroupFunction}
\alias{getFitUnitGroupFunction.MultiArrayUnitModel}
\alias{MultiArrayUnitModel.getFitUnitGroupFunction}
\alias{getFitUnitGroupFunction,MultiArrayUnitModel-method}

\title{Static method to get the low-level function that fits the PLM}

\description{
 Static method to get the low-level function that fits the PLM.
 Any subclass model must provide this method, which should return
 a \code{\link[base]{function}} that accepts an IxK \code{\link[base]{matrix}}.
}

\usage{
## Static method (use this):
## MultiArrayUnitModel$getFitUnitGroupFunction()

## Don't use the below:
\method{getFitUnitGroupFunction}{MultiArrayUnitModel}(...)
}

\arguments{
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{function}}.
}

\seealso{
  For more information see \code{\link{MultiArrayUnitModel}}.
}
\keyword{internal}
\keyword{methods}
