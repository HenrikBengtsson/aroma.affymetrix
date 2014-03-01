%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractProbeSequenceNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AbstractProbeSequenceNormalization}
\docType{class}
\alias{AbstractProbeSequenceNormalization}


\title{The AbstractProbeSequenceNormalization class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class AbstractProbeSequenceNormalization}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform3}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AbstractProbeSequenceNormalization}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{BaseCountNormalization}, \emph{\link[aroma.affymetrix]{BasePositionNormalization}}, \emph{\link[aroma.affymetrix]{LinearModelProbeSequenceNormalization}}, \link[aroma.affymetrix]{MatNormalization}\cr

public abstract static class \bold{AbstractProbeSequenceNormalization}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelTransform3}}\cr



 This abstract class represents a normalization method that corrects for
 systematic effects in the probe intensities due to differences in
 probe sequences.
}

\usage{
AbstractProbeSequenceNormalization(..., target=NULL)
}

\arguments{
  \item{...}{Arguments passed to the constructor of
    \code{\link{ProbeLevelTransform3}}.}
  \item{target}{A \code{\link[base]{character}} string specifying type of "target" used.
    If \code{"zero"}, all arrays are normalized to have no effects.
    If \code{\link[base]{NULL}}, all arrays a normalized to have the same effect as
    the average array has.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getTargetFile} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:process.AbstractProbeSequenceNormalization]{process}} \tab Normalizes the data set.\cr
}


 \bold{Methods inherited from ProbeLevelTransform3}:\cr
getAsteriskTags, getCellsTo, getCellsToFit, getCellsToUpdate, getParameters, getUnitsTo, getUnitsToFit, getUnitsToUpdate, writeSignals

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

\section{Requirements}{
  This class requires that an \code{\link[aroma.core]{AromaCellSequenceFile}} is
  available for the chip type.
}

\author{Henrik Bengtsson}
\keyword{classes}
