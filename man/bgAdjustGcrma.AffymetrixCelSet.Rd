%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.BG.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{bgAdjustGcrma.AffymetrixCelSet}
\alias{bgAdjustGcrma.AffymetrixCelSet}
\alias{AffymetrixCelSet.bgAdjustGcrma}
\alias{bgAdjustGcrma,AffymetrixCelSet-method}

\title{Applies probe sequence based background correction to a set of
CEL files}

\description{
 Applies probe sequence based background correction to a set of
CEL files.

 Adapted from \code{\link[gcrma]{bg.adjust.gcrma}} in the \pkg{gcrma} package.
}

\usage{
\method{bgAdjustGcrma}{AffymetrixCelSet}(this, path, affinities=NULL, type="fullmodel", indicesNegativeControl=NULL,
  opticalAdjust=TRUE, gsbAdjust=TRUE, k=6 * fast + 0.5 * (1 - fast), rho=0.7,
  stretch=1.15 * fast + (1 - fast), fast=TRUE, overwrite=FALSE, skip=!overwrite, ...,
  verbose=FALSE, .deprecated=TRUE)
}

\arguments{
  \item{path}{The path where to save the adjusted data files.}
  \item{name}{Name of the set containing the background corrected files.}
  \item{type}{The type of background correction.  Currently accepted types
      are "fullmodel" (the default, uses MMs) and "affinities" (uses
      probe sequence only).}
  \item{indicesNegativeControl}{Locations of any negative control
      probes (e.g., the anti-genomic controls on the human exon array).
      If \code{\link[base]{NULL}} and type=="affinities", MMs are used as the negative
      controls.}
  \item{opticalAdjust}{If \code{\link[base:logical]{TRUE}}, apply correction for optical effect,
      as in \code{\link[gcrma]{bg.adjust.optical}}.}
  \item{gsbAdjust}{Should we adjust for specific binding (defaults to
       \code{\link[base:logical]{TRUE}})?}
  \item{k}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
  \item{rho}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
  \item{stretch}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
  \item{fast}{If \code{\link[base:logical]{TRUE}}, an ad hoc transformation of the PM is performed
      (\code{gcrma::gcrma.bg.transformation.fast}).}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already adjusted arrays are overwritten,
    unless skipped, otherwise an error is thrown.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the array is not normalized if it already exists.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{.deprecated}{Internal argument.}
}

\value{
 Returns the background adjusted \code{\link{AffymetrixCelFile}} object.
}

\author{Ken Simpson}

\seealso{
 \code{\link[gcrma]{bg.adjust.gcrma}}
 For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
