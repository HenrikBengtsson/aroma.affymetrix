%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AffymetrixCelFile$clearData}
\alias{AffymetrixCelFile$clearData}
\alias{clearData.AffymetrixCelFile}
\alias{AffymetrixCelFile.clearData}
\alias{clearData,AffymetrixCelFile-method}

\title{Clears all or a subset of the fields in a CEL file}

\description{
 Clears all or a subset of the fields in a CEL file.
}

\usage{
## Static method (use this):
## AffymetrixCelFile$clearData(fields=c("intensities", "stdvs", "pixels"), value=0, ...,
##   .forSure=FALSE, verbose=TRUE)

## Don't use the below:
\method{clearData}{AffymetrixCelFile}(this, fields=c("intensities", "stdvs", "pixels"), value=0, ..., .forSure=FALSE,
  verbose=TRUE)
}

\arguments{
  \item{fields}{A \code{\link[base]{character}} \code{\link[base]{vector}} of fields to be cleared.}
  \item{value}{A \code{\link[base]{numeric}} value to be written over the data.}
  \item{...}{Not used.}
  \item{.forSure}{If not \code{\link[base:logical]{TRUE}}, an exception is thrown asking if the
     method was called by mistake.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns (invisibly) the names of the fields cleared.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
  Internally, \code{\link[affxparser]{updateCel}} is used.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
