%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectFile.fromDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{ChipEffectFile$fromDataFile}
\alias{ChipEffectFile$fromDataFile}
\alias{fromDataFile.ChipEffectFile}
\alias{ChipEffectFile.fromDataFile}
\alias{fromDataFile,ChipEffectFile-method}

\title{Retrieve an existing CEL file, or create from CDF template if missing}

\description{
 Retrieve an existing CEL file, or create from CDF template if missing.
}

\usage{
## Static method (use this):
## ChipEffectFile$fromDataFile(df=NULL, filename=sprintf("\%s,chipEffects.CEL",
##   getFullName(df)), path, name=getName(df), cdf=NULL, ..., verbose=FALSE)

## Don't use the below:
\method{fromDataFile}{ChipEffectFile}(static, df=NULL, filename=sprintf("\%s,chipEffects.CEL", getFullName(df)),
  path, name=getName(df), cdf=NULL, ..., verbose=FALSE)
}

\arguments{
  \item{static}{}
  \item{df}{}
  \item{filename}{The filename of the CEL file.}
  \item{path}{The path to the directory where to find/create the CEL file.}
  \item{name}{The name of the array to be stored in the CEL header.}
  \item{cdf}{The template \code{\link{AffymetrixCdfFile}} used for creating
             a CEL file from scratch.}
  \item{...}{Passed to \code{\link[affxparser]{createCel}}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link{ChipEffectFile}}.
}

\author{Henrik Bengtsson}

\seealso{
  \code{allocateFromCdf()} of \code{\link{AffymetrixCelFile}}.
  For more information see \code{\link{ChipEffectFile}}.
}
\keyword{internal}
\keyword{methods}
