%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ReseqCrosstalkCalibration.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ReseqCrosstalkCalibration}
\docType{class}
\alias{ReseqCrosstalkCalibration}


\title{The ReseqCrosstalkCalibration class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class ReseqCrosstalkCalibration}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{ReseqCrosstalkCalibration}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ReseqCrosstalkCalibration}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr



 This class represents a calibration function that transforms the
 probe-level signals such that the signals from the four nucleotides
 (A, C, G, T) are orthogonal.
}

\usage{
ReseqCrosstalkCalibration(dataSet=NULL, ..., targetAvg=2200, subsetToAvg=NULL,
  mergeGroups=FALSE, flavor=c("sfit", "expectile"), alpha=c(0.1, 0.075, 0.05, 0.03,
  0.01), q=2, Q=98)
}

\arguments{
  \item{dataSet}{An \code{\link{AffymetrixCelSet}}.}
  \item{...}{Arguments passed to the constructor of
    \code{\link{ProbeLevelTransform}}.}
  \item{targetAvg}{The signal(s) that the average of the sum of the
    probe quartets should have after calibration.}
  \item{subsetToAvg}{The indices of the cells (taken as the intersect of
    existing indices) used to calculate average in order to rescale to
    the target average. If \code{\link[base]{NULL}}, all probes are considered.}
  \item{mergeGroups}{A \code{\link[base]{logical}} ...}
  \item{flavor}{A \code{\link[base]{character}} string specifying what algorithm is used
    to fit the crosstalk calibration.}
  \item{alpha, q, Q}{Additional arguments passed to
   \code{fitMultiDimensionalCone()}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.ReseqCrosstalkCalibration]{process}} \tab Calibrates the data set.\cr
}


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

\author{Henrik Bengtsson}
\keyword{classes}
