%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectFile.fromDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fromDataFile.ChipEffectFile}
\alias{fromDataFile.ChipEffectFile}
\alias{ChipEffectFile.fromDataFile}
\alias{fromDataFile.ChipEffectFile}
\alias{fromDataFile,ChipEffectFile-method}


\title{Retrive an existing CEL file, or create from CDF template if missing}

\description{
 Retrive an existing CEL file, or create from CDF template if missing.
}

\synopsis{fromDataFile.ChipEffectFile(static, df=NULL, filename=sprintf("\%s,chipEffects.CEL", getFullName(df)), path, name=getName(df), cdf=NULL, ..., verbose=FALSE)}
\usage{ChipEffectFile$fromDataFile(df=NULL, filename=sprintf("\%s,chipEffects.CEL", getFullName(df)), path, name=getName(df), cdf=NULL, ..., verbose=FALSE)}

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

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{allocateFromCdf()} of \code{\link{AffymetrixCelFile}}.
  For more information see \code{\link{ChipEffectFile}}.
}
\keyword{internal}
\keyword{methods}
