%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCellIndices.AffymetrixCdfFile}
\alias{getCellIndices.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.getCellIndices}
\alias{getCellIndices.AffymetrixCdfFile}
\alias{getCellIndices,AffymetrixCdfFile-method}


\title{Gets the cell indices unit by unit}

\description{
 Gets the cell indices unit by unit of all or a subset of the units.
}

\usage{\method{getCellIndices}{AffymetrixCdfFile}(this, units=NULL, ..., useNames=TRUE, unlist=FALSE, force=FALSE, cache=TRUE, verbose=FALSE)}

\arguments{
  \item{units}{The units of interest. If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Additional arguments passed to
     \code{\link[affxparser]{readCdfCellIndices}}.}
  \item{useNames}{If \code{\link[base:logical]{TRUE}}, element names are returned, otherwise not.}
  \item{unlist}{If \code{\link[base:logical]{TRUE}}, the unlisted result is returned. Using this
     argument is more memory efficient that calling \code{\link[base]{unlist}}()
     afterwards.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, cached values are ignored.}
  \item{cache}{If \code{\link[base:logical]{TRUE}}, results are cached, if not too large.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the \code{\link[base]{list}} structure returned by
 \code{\link[affxparser]{readCdfCellIndices}}.
}

\details{
 Note, that it is much more memory efficient to do
 \code{getCellIndices(cdf, useNames=FALSE, unlist=TRUE)}
 compare with \code{unlist(getCellIndices(cdf), use.names=FALSE)}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  See \code{\link[aroma.affymetrix:setRestructor.AffymetrixCdfFile]{*setRestructor}()} to set a default re-constructor for
  the returned CDF structure.
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
