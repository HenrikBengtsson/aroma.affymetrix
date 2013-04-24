%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaChipTypeAnnotationFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findByChipType.AromaChipTypeAnnotationFile}
\alias{findByChipType.AromaChipTypeAnnotationFile}
\alias{AromaChipTypeAnnotationFile.findByChipType}
\alias{findByChipType,AromaChipTypeAnnotationFile-method}

\title{Locates an annotation file by its chip type}

\description{
 Locates an annotation file by its chip type.
}

\usage{
\method{findByChipType}{AromaChipTypeAnnotationFile}(...)
}

\arguments{
 \item{chipType}{A \code{\link[base]{character}} string.}
 \item{tags}{An optional \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
 \item{...}{Additional arguments.}
}

\value{
 Returns the pathname (as a \code{\link[base]{character}} string) to the first annotation
 chip type file found.  If no one was found, \code{\link[base]{NULL}} is returned.
}

\seealso{
  \code{\link[aroma.affymetrix:byChipType.AromaChipTypeAnnotationFile]{*byChipType}()}.
  For more information see \code{\link{AromaChipTypeAnnotationFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
