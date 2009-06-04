%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.BG.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{bgAdjustGcrma.AffymetrixCelFile}
\alias{bgAdjustGcrma.AffymetrixCelFile}
\alias{AffymetrixCelFile.bgAdjustGcrma}
\alias{bgAdjustGcrma.AffymetrixCelFile}
\alias{bgAdjustGcrma,AffymetrixCelFile-method}


\title{Applies probe sequence based background correction to a CEL file}

\description{
 Applies probe sequence based background correction to a CEL file.

 Adapted from \code{\link[gcrma]{bg.adjust.gcrma}} in the \pkg{gcrma} package.
}

\usage{\method{bgAdjustGcrma}{AffymetrixCelFile}(this, path=NULL, type="fullmodel", indicesNegativeControl=NULL, affinities=NULL, gsbAdjust=TRUE, parametersGsb=NULL, k=6 * fast + 0.5 * (1 - fast), rho=0.7, stretch=1.15 * fast + 1 * (1 - fast), fast=TRUE, overwrite=FALSE, skip=!overwrite, ..., verbose=FALSE, .deprecated=TRUE)}

\arguments{
  \item{path}{The path where to save the adjusted data files.}
  \item{type}{The type of background correction.  Currently accepted types
      are "fullmodel" (the default, uses MMs) and "affinities" (uses
      probe sequence only).}
  \item{indicesNegativeControl}{Locations of any negative control
      probes (e.g., the anti-genomic controls on the human exon array).
      If \code{\link[base]{NULL}} and type=="affinities", MMs are used as the negative
      controls.}
  \item{affinities}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of probe affinities, usually as
      calculated by \code{computeAffinities()} of the
      \code{\link{AffymetrixCdfFile}} class.}
  \item{gsbAdjust}{Should we adjust for specific binding (defaults to
       \code{\link[base:logical]{TRUE}})?}
  \item{parametersGsb}{Specific binding parameters as estimated by
      \code{calculateParametersGsb()} for the \code{\link{AffymetrixCelSet}}
      class.}
  \item{k}{Tuning parameter passed to \code{bg.adjust.gcrma}.}
  \item{rho}{Tuning parameter passed to \code{bg.adjust.gcrma}.}
  \item{stretch}{Tuning parameter passed to \code{bg.adjust.gcrma}.}
  \item{fast}{If \code{\link[base:logical]{TRUE}}, an ad hoc transformation of the PM is performed
      (\code{gcrma.bg.transformation.fast()}).}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already adjusted arrays are overwritten,
    unless skipped, otherwise an error is thrown.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the array is not normalized if it already exists.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{.deprecated}{Internal argument.}
}

\value{
 Returns the background adjusted \code{\link{AffymetrixCelFile}} object.
}

\author{
  Ken Simpson (ksimpson[at]wehi.edu.au) and Mark Robinson.
}

\seealso{
 \code{\link[gcrma]{bg.adjust.gcrma}}
 For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
