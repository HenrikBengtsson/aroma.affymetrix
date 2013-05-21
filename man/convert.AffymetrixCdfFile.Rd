%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{convert.AffymetrixCdfFile}
\alias{convert.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.convert}
\alias{convert,AffymetrixCdfFile-method}

\title{Converts a CDF into the same CDF but with another format}

\description{
 Converts a CDF into the same CDF but with another format.
}

\usage{
\method{convert}{AffymetrixCdfFile}(this, chipType=getChipType(this), suffix=NULL, sep="-", path="cdf", ...,
  verbose=FALSE)
}

\arguments{
  \item{chipType}{The chip type of the new CDF.}
  \item{suffix}{A suffix added to the chip type of the new CDF.}
  \item{sep}{A string separating the chip type and the suffix string.}
  \item{path}{The path where to store the new CDF file.}
  \item{...}{Additional arguments passed to \code{\link[affxparser]{convertCdf}}.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the new CDF as an \code{\link{AffymetrixCdfFile}} object.
}

\seealso{
  To compare two CDFs, use \code{equals()}.
  Internally \code{\link[affxparser]{convertCdf}} is used.
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
