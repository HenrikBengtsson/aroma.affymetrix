%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChromosomalModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getTableOfArrays.ChromosomalModel}
\alias{getTableOfArrays.ChromosomalModel}
\alias{ChromosomalModel.getTableOfArrays}
\alias{getTableOfArrays.ChromosomalModel}
\alias{getTableOfArrays,ChromosomalModel-method}


\title{Gets a table of arrays}

\description{
 Gets a table of arrays showing their availability across chip types.
}

\usage{\method{getTableOfArrays}{ChromosomalModel}(this, ...)}

\arguments{
  \item{...}{Not used.}
}

\value{
 Returns a \eqn{NxK} \code{\link[base]{matrix}} of \code{\link[base]{integer}}s where \eqn{N} is the total number
 of arrays and \eqn{K} is the number of chip types in the model.  The row
 names are the names of the arrays, and the column names are the chip types.
 If data is available for array \eqn{n} and chip type \eqn{k}, cell
 \eqn{(n,k)} has value \eqn{n}, otherwise \code{\link[base]{NA}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ChromosomalModel}}.
}
\keyword{internal}
\keyword{methods}
