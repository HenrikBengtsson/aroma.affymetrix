%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ParameterCelSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{extractDataFrame.ParameterCelSet}
\alias{extractDataFrame.ParameterCelSet}
\alias{ParameterCelSet.extractDataFrame}
\alias{extractDataFrame.ParameterCelSet}
\alias{extractDataFrame,ParameterCelSet-method}


\title{Extract data as a data.frame for a set of arrays}

\description{
 Extract data as a data.frame for a set of arrays.
}

\usage{\method{extractDataFrame}{ParameterCelSet}(this, addNames=FALSE, addUgcMap=TRUE, ..., drop=FALSE, verbose=FALSE)}

\arguments{
  \item{...}{Arguments passed to \code{\link[aroma.affymetrix:extractMatrix.ParameterCelSet]{*extractMatrix}()}.}
  \item{addNames}{If \code{\link[base:logical]{TRUE}}, the first two columns contain the
    unit names and the group names according the the CDF, otherwise
    those two columns are not included.}
  \item{addUgcMap}{If \code{\link[base:logical]{TRUE}}, the columns following the unit and
    group names contains the (unit, group, cell) index map.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a Jx(2+3+K) \code{\link[base]{data.frame}} where J is the number of units,
 and K is the number of arrays.
 The first two columns, if \code{addNames=TRUE}, contains the
 unit names and the group names.
 The next three columns contains the (unit, group, cell) index map.
 The last K columns named by the arrays contain the data for the K arrays.
 No names are set for the rows.
 The rows are ordered according to \code{units} arguments.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:extractMatrix.ParameterCelSet]{*extractMatrix}()}.
  For more information see \code{\link{ParameterCelSet}}.
}
\keyword{internal}
\keyword{methods}
