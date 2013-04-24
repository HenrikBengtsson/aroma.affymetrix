%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectGroupMerge.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ChipEffectGroupMerge}
\docType{class}
\alias{ChipEffectGroupMerge}


\title{The ChipEffectGroupMerge class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class ChipEffectGroupMerge}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectTransform}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{ChipEffectGroupMerge}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{SnpChipEffectGroupMerge}\cr

public static class \bold{ChipEffectGroupMerge}\cr
extends \emph{\link[aroma.affymetrix]{ChipEffectTransform}}\cr



 This class represents a method that merges chip effects across groups
 unit by unit.
}

\usage{
ChipEffectGroupMerge(dataSet=NULL, fcn=NULL, ...)
}

\arguments{
  \item{dataSet}{A \code{\link{ChipEffectSet}}.}
  \item{fcn}{A \code{\link[base]{function}} used to "merge" groups.}
  \item{...}{Additional arguments passed to the constructor of
    \code{\link{ChipEffectTransform}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getMergeFunction} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:process.ChipEffectGroupMerge]{process}} \tab Normalizes the data set.\cr
 \tab \code{setMergeFunction} \tab  -\cr
}


 \bold{Methods inherited from ChipEffectTransform}:\cr
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
\keyword{internal}
