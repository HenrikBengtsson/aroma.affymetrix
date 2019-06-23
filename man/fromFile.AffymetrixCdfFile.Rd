%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AffymetrixCdfFile$fromFile}
\alias{AffymetrixCdfFile$fromFile}
\alias{fromFile.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.fromFile}
\alias{fromFile,AffymetrixCdfFile-method}

\title{Defines an AffymetrixCdfFile object from a CDF file}

\description{
 Defines an AffymetrixCdfFile object from a CDF file.
}

\usage{
## Static method (use this):
## AffymetrixCdfFile$fromFile(filename, path=NULL, ...)

## Don't use the below:
\method{fromFile}{AffymetrixCdfFile}(static, filename, path=NULL, ...)
}

\arguments{
 \item{filename}{The filename of to the file.}
 \item{path}{The path to the file.}
 \item{...}{Not used.}
}

\value{
 Returns an instance of \code{\link{AffymetrixCdfFile}} or its subclasses.
 If the file is not found or if it is of the wrong file format, an
 error is thrown.
}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
