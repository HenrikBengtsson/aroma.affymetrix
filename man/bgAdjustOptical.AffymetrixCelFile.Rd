%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.BG.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{bgAdjustOptical.AffymetrixCelFile}
\alias{bgAdjustOptical.AffymetrixCelFile}
\alias{AffymetrixCelFile.bgAdjustOptical}
\alias{bgAdjustOptical,AffymetrixCelFile-method}

\title{Applies optical background correction to a CEL file}

\description{
 Applies optical background correction to a CEL file.

Adapted from \code{\link[gcrma]{bg.adjust.optical}} in the \pkg{gcrma} package.
}

\usage{
\method{bgAdjustOptical}{AffymetrixCelFile}(this, path, minimum=1, subsetToUpdate=NULL, typesToUpdate=NULL,
  overwrite=FALSE, skip=!overwrite, verbose=FALSE, ..., .deprecated=TRUE)
}

\arguments{
  \item{path}{The location to save the adjusted data files.}
  \item{minimum}{The minimum adjusted intensity.  Defaults to 1.}
  \item{subsetToUpdate}{The indices of the probes to be updated.
    If \code{\link[base]{NULL}}, all are updated.}
  \item{typesToUpdate}{Types of probes to be updated.  For more details,
    see argument \code{types} of \code{identifyCells()} for the
    \code{\link{AffymetrixCdfFile}} class.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already adjusted arrays are overwritten,
    unless skipped, otherwise an error is thrown.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the array is not normalized if it already exists.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{.deprecated}{Internal argument.}
}

\value{
 Returns the background adjusted \code{\link{AffymetrixCelFile}} object.
}

\author{Ken Simpson, Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
