%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromFile.AffymetrixCdfFile}
\alias{fromFile.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.fromFile}
\alias{fromFile.AffymetrixCdfFile}
\alias{fromFile,AffymetrixCdfFile-method}


\title{Defines an AffymetrixCdfFile object from a CDF file}

\description{
 Defines an AffymetrixCdfFile object from a CDF file.
}

\synopsis{fromFile.AffymetrixCdfFile(static, filename, path=NULL, ...)}
\usage{AffymetrixCdfFile$fromFile(filename, path=NULL, ...)}

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

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:byChipType.AffymetrixCdfFile]{*byChipType}()}.
  For more information see \code{\link{AffymetrixCdfFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
