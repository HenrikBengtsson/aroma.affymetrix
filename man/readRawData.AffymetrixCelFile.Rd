%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readRawData.AffymetrixCelFile}
\alias{readRawData.AffymetrixCelFile}
\alias{AffymetrixCelFile.readRawData}
\alias{readRawData,AffymetrixCelFile-method}
\alias{AffymetrixCelFile.getData}
\alias{getData.AffymetrixCelFile}
\alias{getData,AffymetrixCelFile-method}


\title{Gets all or a subset of the fields in a CEL file}

\description{
 Gets all or a subset of the fields in a CEL file for all or a subset of the cells.
}

\usage{
\method{readRawData}{AffymetrixCelFile}(this, indices=NULL, fields=c("xy", "intensities", "stdvs", "pixels"), ...,
  drop=FALSE, verbose=FALSE)
}

\arguments{
  \item{indices}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of cell indices.  If \code{\link[base]{NULL}}, all cells
    are considered.}
  \item{fields}{Names of fields to be retrieved.}
  \item{...}{Additional arguments passed to \code{\link[affxparser]{readCel}}.}
  \item{drop}{If \code{\link[base:logical]{TRUE}} and a single field is returned, then data is
    returned as a \code{\link[base]{vector}}, otherwise as a \code{\link[base]{data.frame}}.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{data.frame}} of the fields requested (unless dimension dropped).
}

\section{Caching}{
  Neither in-memory nor on-file caching is done by this method.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
