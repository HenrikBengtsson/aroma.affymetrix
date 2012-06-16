%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixPgfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findByChipType.AffymetrixPgfFile}
\alias{findByChipType.AffymetrixPgfFile}
\alias{AffymetrixPgfFile.findByChipType}
\alias{findByChipType.AffymetrixPgfFile}
\alias{findByChipType,AffymetrixPgfFile-method}


\title{Locates a PGF file from its chip type}

\description{
 Locates a PGF file from its chip type.
}

\synopsis{findByChipType.AffymetrixPgfFile(static, chipType, tags=NULL, pattern=NULL, ...)}
\usage{AffymetrixPgfFile$findByChipType(chipType, tags=NULL, pattern=NULL, ...)}

\arguments{
 \item{chipType}{A \code{\link[base]{character}} string.}
 \item{tags}{An optional \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
 \item{pattern}{An optional \code{\link[base]{character}} string.}
 \item{...}{Not used.}
}

\value{
 Returns a pathname as a \code{\link[base]{character}} string to the first PGF file found.
 If non PGF with requested chip type was found, \code{\link[base]{NULL}} is returned.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixPgfFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
