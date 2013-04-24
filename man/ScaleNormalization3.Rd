%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ScaleNormalization3.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ScaleNormalization3}
\docType{class}
\alias{ScaleNormalization3}


\title{The ScaleNormalization3 class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class ScaleNormalization3}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{ScaleNormalization3}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ScaleNormalization3}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelTransform3}}\cr



 This class represents a normalization function that transforms the
 probe-level signals towards the same scale.
}

\usage{
ScaleNormalization3(..., targetAvg=4400)
}

\arguments{
  \item{...}{Arguments passed to the constructor of
    \code{\link{ProbeLevelTransform3}}.}
  \item{targetAvg}{A \code{\link[base]{numeric}} value.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.ScaleNormalization3]{process}} \tab Normalizes the data set.\cr
}


 \bold{Methods inherited from ProbeLevelTransform3}:\cr
getAsteriskTags, getCellsTo, getCellsToFit, getCellsToUpdate, getParameters, getUnitsTo, getUnitsToFit, getUnitsToUpdate, writeSignals

 \bold{Methods inherited from ProbeLevelTransform}:\cr
getRootPath

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
as.character, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
