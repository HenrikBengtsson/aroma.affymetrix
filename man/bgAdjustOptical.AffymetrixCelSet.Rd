%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.BG.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{bgAdjustOptical.AffymetrixCelSet}
\alias{bgAdjustOptical.AffymetrixCelSet}
\alias{AffymetrixCelSet.bgAdjustOptical}
\alias{bgAdjustOptical.AffymetrixCelSet}
\alias{bgAdjustOptical,AffymetrixCelSet-method}


\title{Applies optical background correction to a set of CEL files}

\description{
 Applies optical background correction to a set of CEL files.

 Adapted from \code{\link[gcrma]{bg.adjust.optical}} in the \pkg{gcrma} package.
}

\usage{\method{bgAdjustOptical}{AffymetrixCelSet}(this, path=NULL, name="bgOptical", subsetToUpdate=NULL, typesToUpdate=NULL, minimum=1, overwrite=FALSE, skip=!overwrite, ..., verbose=FALSE, .deprecated=TRUE)}

\arguments{
  \item{path}{The location to save the adjusted data files.}
  \item{minimum}{The minimum adjusted intensity.  Defaults to 1.}
  \item{subsetToUpdate}{The indices of the probes to be updated.
    If \code{\link[base]{NULL}}, all are updated.}
  \item{typesToUpdate}{Types of probes to be updated.  For more details,
    see argument \code{types} of \code{identifyCells()} for the
    \code{\link{AffymetrixCdfFile}} class.}
  \item{...}{Not used.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already adjusted arrays are overwritten,
    unless skipped, otherwise an error is thrown.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the array is not normalized if it already exists.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{.deprecated}{Internal argument.}
}

\value{
 Returns the background adjusted \code{\link{AffymetrixCelSet}} object.
}

\author{
  Ken Simpson (ksimpson[at]wehi.edu.au).
}
\keyword{internal}
\keyword{methods}
