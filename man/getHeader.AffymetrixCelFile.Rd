%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getHeader.AffymetrixCelFile}
\alias{getHeader.AffymetrixCelFile}
\alias{AffymetrixCelFile.getHeader}
\alias{getHeader,AffymetrixCelFile-method}

\title{Gets the header of the CEL file}

\description{
 Gets the header of the CEL file.
}

\usage{
\method{getHeader}{AffymetrixCelFile}(this, ...)
}

\arguments{
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{list}} structure as returned by \code{\link[affxparser]{readCelHeader}}.
 If file does not exists, then \code{\link[base]{NULL}} is returned.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
