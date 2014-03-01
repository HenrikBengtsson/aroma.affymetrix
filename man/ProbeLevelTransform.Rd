%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeLevelTransform.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ProbeLevelTransform}
\docType{class}
\alias{ProbeLevelTransform}


\title{The ProbeLevelTransform class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class ProbeLevelTransform}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\emph{\code{ProbeLevelTransform}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AbstractProbeSequenceNormalization}}, \link[aroma.affymetrix]{AllelicCrosstalkCalibration}, \emph{\link[aroma.affymetrix]{BackgroundCorrection}}, \link[aroma.affymetrix]{BaseCountNormalization}, \emph{\link[aroma.affymetrix]{BasePositionNormalization}}, \link[aroma.affymetrix]{DChipQuantileNormalization}, \link[aroma.affymetrix]{GcRmaBackgroundCorrection}, \link[aroma.affymetrix]{LimmaBackgroundCorrection}, \emph{\link[aroma.affymetrix]{LinearModelProbeSequenceNormalization}}, \link[aroma.affymetrix]{MatNormalization}, \link[aroma.affymetrix]{MatSmoothing}, \link[aroma.affymetrix]{NormExpBackgroundCorrection}, \link[aroma.affymetrix]{OpticalBackgroundCorrection}, \emph{\link[aroma.affymetrix]{ProbeLevelTransform3}}, \link[aroma.affymetrix]{QuantileNormalization}, \link[aroma.affymetrix]{ReseqCrosstalkCalibration}, \link[aroma.affymetrix]{RmaBackgroundCorrection}, \link[aroma.affymetrix]{ScaleNormalization}, \link[aroma.affymetrix]{ScaleNormalization3}, \link[aroma.affymetrix]{SpatialRowColumnNormalization}, \link[aroma.affymetrix]{UnitTypeScaleNormalization}\cr

public abstract static class \bold{ProbeLevelTransform}\cr
extends \emph{\link[aroma.affymetrix]{Transform}}\cr



 This abstract class represents a transformation methods that transforms
 probe-level signals, typically intensities.
}

\usage{
ProbeLevelTransform(...)
}

\arguments{
  \item{...}{Arguments passed to the constructor of \code{\link{Transform}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
as.character, findFilesTodo, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, asThis


}

\details{
  Subclasses must implement the \code{process()} method.
}

\author{Henrik Bengtsson}
\keyword{classes}
