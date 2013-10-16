%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeLevelTransform3.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ProbeLevelTransform3}
\docType{class}
\alias{ProbeLevelTransform3}


\title{The ProbeLevelTransform3 class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class ProbeLevelTransform3}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{ProbeLevelTransform3}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AbstractProbeSequenceNormalization}}, \link[aroma.affymetrix]{BaseCountNormalization}, \emph{\link[aroma.affymetrix]{BasePositionNormalization}}, \emph{\link[aroma.affymetrix]{LinearModelProbeSequenceNormalization}}, \link[aroma.affymetrix]{MatNormalization}, \link[aroma.affymetrix]{ScaleNormalization3}, \link[aroma.affymetrix]{UnitTypeScaleNormalization}\cr

public abstract static class \bold{ProbeLevelTransform3}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr



 This abstract class is specialized from \code{\link{ProbeLevelTransform}} and
 provides methods to identify subsets and types of probes that are used
 for fitting and/or updating the signals.
}

\usage{
ProbeLevelTransform3(dataSet=NULL, ..., unitsToFit="-XY", typesToFit=typesToUpdate,
  unitsToUpdate=NULL, typesToUpdate="pm", shift=0)
}

\arguments{
  \item{dataSet}{A \code{\link{AffymetrixCelSet}}.}
  \item{...}{Arguments passed to the constructor of
    \code{\link{ProbeLevelTransform}}.}
  \item{unitsToFit}{The units from which the normalization curve should
    be estimated.  If \code{\link[base]{NULL}}, all are considered.}
  \item{typesToFit}{Types of probes to be used when fitting the model.}
  \item{unitsToUpdate}{The units to be updated.
    If \code{\link[base]{NULL}}, all are considered.}
  \item{typesToUpdate}{Types of probes to be updated.}
  \item{shift}{An optional amount to shift data before fitting and updating.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from ProbeLevelTransform}:\cr
getRootPath

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
as.character, findFilesTodo, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
