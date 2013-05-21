%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.BG.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{bgAdjustRma.AffymetrixCelFile}
\alias{bgAdjustRma.AffymetrixCelFile}
\alias{AffymetrixCelFile.bgAdjustRma}
\alias{bgAdjustRma,AffymetrixCelFile-method}

\title{Applies normExp background correction to a CEL file}

\description{
 Applies normExp background correction to a CEL file.

 Calls \code{\link[affy]{bg.adjust}} from the \pkg{affy} package.
}

\usage{
\method{bgAdjustRma}{AffymetrixCelFile}(this, path, pmonly=TRUE, addJitter=FALSE, jitterSd=0.2, overwrite=FALSE,
  skip=!overwrite, ..., verbose=FALSE, .deprecated=TRUE)
}

\arguments{
  \item{path}{The path where to save the adjusted data files.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already adjusted arrays are overwritten,
    unless skipped, otherwise an error is thrown.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the array is not normalized if it already exists.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{...}{Not used.}
  \item{.deprecated}{Internal argument.}
}

\value{
 Returns the background adjusted \code{\link{AffymetrixCelFile}} object.
}

\author{Ken Simpson, Henrik Bengtsson}

\seealso{
 \code{\link[affy]{bg.adjust}}
 For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
