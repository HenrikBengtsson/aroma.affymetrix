%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultiArrayUnitModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getFitUnitGroupFunction.MultiArrayUnitModel}
\alias{getFitUnitGroupFunction.MultiArrayUnitModel}
\alias{MultiArrayUnitModel.getFitUnitGroupFunction}
\alias{getFitUnitGroupFunction.MultiArrayUnitModel}
\alias{getFitUnitGroupFunction,MultiArrayUnitModel-method}


\title{Static method to get the low-level function that fits the PLM}

\description{
 Static method to get the low-level function that fits the PLM.
 Any subclass model must provide this method, which should return
 a \code{\link[base]{function}} that accepts an IxK \code{\link[base]{matrix}}.
}

\synopsis{getFitUnitGroupFunction.MultiArrayUnitModel(...)}
\usage{MultiArrayUnitModel$getFitUnitGroupFunction()}

\arguments{
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{function}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MultiArrayUnitModel}}.
}
\keyword{internal}
\keyword{methods}
