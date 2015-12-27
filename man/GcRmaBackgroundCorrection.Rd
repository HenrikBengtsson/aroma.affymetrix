%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GcRmaBackgroundCorrection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{GcRmaBackgroundCorrection}
\docType{class}
\alias{GcRmaBackgroundCorrection}


\title{The GcRmaBackgroundCorrection class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class GcRmaBackgroundCorrection}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{BackgroundCorrection}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{GcRmaBackgroundCorrection}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{GcRmaBackgroundCorrection}\cr
extends \emph{\link[aroma.affymetrix]{BackgroundCorrection}}\cr



 This class represents the GCRMA background adjustment function.

}

\usage{
GcRmaBackgroundCorrection(..., indicesNegativeControl=NULL, affinities=NULL,
  type=c("fullmodel", "affinities"), opticalAdjust=TRUE, gsbAdjust=TRUE,
  gsbParameters=NULL, seed=NULL)
}

\arguments{
  \item{...}{Arguments passed to the constructor of
     \code{\link{ProbeLevelTransform}}.}
  \item{indicesNegativeControl}{Locations of any negative control
      probes (e.g., the anti-genomic controls on the human exon array).
      If \code{\link[base]{NULL}} and \code{type == "affinities"}, then all non-PM probes
      are used as the negative controls.}
  \item{affinities}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of probe affinities, usually as
      calculated by \code{computeAffinities()} of the
      \code{\link{AffymetrixCdfFile}} class.}
  \item{type}{Type (flavor) of background correction, which can
      be either \code{"fullmodel"} (uses MMs; requires that the chip type
      has PM/MM pairs) or \code{"affinities"} (uses probe sequence only).}
  \item{gsbAdjust}{If \code{\link[base:logical]{TRUE}}, adjustment for specific binding is done,
      otherwise not.}
  \item{opticalAdjust}{If \code{\link[base:logical]{TRUE}}, correction for optical effect is done
      first, utilizing \code{\link{OpticalBackgroundCorrection}}.}
  \item{gsbParameters}{Additional argument passed to the internal
      \code{bgAdjustGcrma()} method.}
  \item{seed}{An (optional) \code{\link[base]{integer}} specifying a temporary random seed
    to be used during processing.  The random seed is set to its original
    state when done.  If \code{\link[base]{NULL}}, it is not set.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from BackgroundCorrection}:\cr
getParameters, process

 \bold{Methods inherited from ProbeLevelTransform}:\cr
getRootPath

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
as.character, findFilesTodo, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, asThis


}

\references{
 [1] Z. Wu, R. Irizarry, R. Gentleman, F.M. Murillo & F. Spencer.
     \emph{A Model Based Background Adjustment for Oligonucleotide
     Expression Arrays}, JASA, 2004.\cr
}

\author{Ken Simpson, Henrik Bengtsson}
\keyword{classes}
