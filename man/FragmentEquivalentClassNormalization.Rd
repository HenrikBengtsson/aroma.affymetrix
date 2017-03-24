%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FragmentEquivalentClassNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FragmentEquivalentClassNormalization}
\docType{class}
\alias{FragmentEquivalentClassNormalization}


\title{The FragmentEquivalentClassNormalization class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class FragmentEquivalentClassNormalization}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{FragmentEquivalentClassNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{FragmentEquivalentClassNormalization}\cr
extends \emph{\link[aroma.affymetrix]{ChipEffectTransform}}\cr



 This class represents a normalization method that corrects for systematic
 effects between loci of different equivalent classes of pairs of
 sequences that are recognized by the restriction enzymes that cut the
 DNA studies.
}

\usage{
FragmentEquivalentClassNormalization(dataSet=NULL, ..., targetAvgs=NULL,
  subsetToFit="-XY")
}

\arguments{
  \item{dataSet}{A \code{\link{CnChipEffectSet}}.}
  \item{...}{Additional arguments passed to the constructor of
    \code{\link{ChipEffectTransform}}.}
  \item{targetAvgs}{An optional list of \code{\link[base]{function}}s.
    For each enzyme there is one target averages to which all arrays
    should be normalized to.}
  \item{subsetToFit}{The units from which the normalization curve should
    be estimated.  If \code{\link[base]{NULL}}, all are considered.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getAromaUfcFile} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from ChipEffectTransform}:\cr
getRootPath

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
as.character, findFilesTodo, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, asThis


}

\section{Requirements}{
  This class requires an UFC (Unit Fragment Class) annotation file.
}

\section{Acknowledgments}{
  The idea of normalization signals stratified on enzyme
  recognition sequences is credited to Jim Veitch and
  Ben Bolstad at Affymetrix Inc. (2008) who have designed
  a similar method for copy number estimation in the
  Affymetrix' Genotype Console v2.
}

\author{Henrik Bengtsson}
\keyword{classes}
