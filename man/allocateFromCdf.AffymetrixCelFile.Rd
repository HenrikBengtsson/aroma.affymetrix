%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.allocateFromCdf.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{allocateFromCdf.AffymetrixCelFile}
\alias{allocateFromCdf.AffymetrixCelFile}
\alias{AffymetrixCelFile.allocateFromCdf}
\alias{allocateFromCdf.AffymetrixCelFile}
\alias{allocateFromCdf,AffymetrixCelFile-method}


\title{Creates an empty CEL file from a template CDF}

\description{
 Creates an empty CEL file from a template CDF.
}

\synopsis{allocateFromCdf.AffymetrixCelFile(static, cdf, name, tags=NULL, path=".", suffix=".CEL", ..., overwrite=FALSE, verbose=FALSE)}
\usage{AffymetrixCelFile$allocateFromCdf(cdf, name, tags=NULL, path=".", suffix=".CEL", ..., overwrite=FALSE, verbose=FALSE)}

\arguments{
  \item{cdf}{A \code{\link{AffymetrixCdfFile}} to be used as a template.}
  \item{name, tags}{The name and the tags of the file created.}
  \item{path}{The directory where the file is created.}
  \item{suffix}{Filename suffix.}
  \item{...}{Arguments passed to \code{\link[affxparser]{createCel}}.}
  \item{overwrite}{If \code{\link[base:logical]{FALSE}} and the file already exists, then an
     error is thrown.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link{AffymetrixCelFile}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
