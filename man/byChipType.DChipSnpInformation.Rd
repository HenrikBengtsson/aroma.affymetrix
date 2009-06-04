%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipSnpInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{byChipType.DChipSnpInformation}
\alias{byChipType.DChipSnpInformation}
\alias{DChipSnpInformation.byChipType}
\alias{byChipType.DChipSnpInformation}
\alias{byChipType,DChipSnpInformation-method}


\title{Defines a DChipSnpInformation object by chip type}

\description{
 Defines a DChipSnpInformation object by chip type.
}

\usage{\method{byChipType}{DChipSnpInformation}(static, chipType, version=NULL, ...)}

\arguments{
 \item{chipType}{A \code{\link[base]{character}} string.}
 \item{rootPath}{A \code{\link[base]{character}} string specifying the root path, i.e.
   the annotation directory.}
 \item{version}{An optional \code{\link[base]{character}} string specifying the version
   string, if more than one version is available.}
 \item{pattern}{An optional filename pattern used to locate the
   dChip genome file.  If \code{\link[base]{NULL}}, a default pattern is used.}
 \item{...}{Not used.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link{DChipSnpInformation}} object.
 If no file was not found, an error is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{DChipSnpInformation}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
