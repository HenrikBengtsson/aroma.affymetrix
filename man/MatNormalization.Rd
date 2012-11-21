%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MatNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{MatNormalization}
\docType{class}
\alias{MatNormalization}


 \title{The MatNormalization class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class MatNormalization}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AbstractProbeSequenceNormalization}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{MatNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{MatNormalization}\cr
extends \emph{\link[aroma.affymetrix]{AbstractProbeSequenceNormalization}}\cr



  This class represents a normalization method that corrects for systematic
  effects in the probe intensities due to differences in the number of
  A, C, G, and T:s and the match scores according to MAT [1].
 }

 \usage{MatNormalization(..., unitsToFit=NULL, model=c("lm"), nbrOfBins=200)}

 \arguments{
   \item{...}{Arguments passed to the constructor of
     \code{\link{AbstractProbeSequenceNormalization}}.}
   \item{unitsToFit}{The units from which the normalization curve should
     be estimated.  If \code{\link[base]{NULL}}, all are considered.}
   \item{model}{A \code{\link[base]{character}} string specifying the model used to fit
     the base-count effects.}
   \item{nbrOfBins}{The number of bins to use for the variance smoothing step.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.MatNormalization]{process}} \tab Normalizes the data set.\cr
}


 \bold{Methods inherited from AbstractProbeSequenceNormalization}:\cr
fitOne, getAromaCellSequenceFile, getParameters, getTargetFile, indexOfMissingSequences, predictOne, process

 \bold{Methods inherited from ProbeLevelTransform3}:\cr
getAsteriskTags, getCellsTo, getCellsToFit, getCellsToUpdate, getParameters, getUnitsTo, getUnitsToFit, getUnitsToUpdate, writeSignals

 \bold{Methods inherited from ProbeLevelTransform}:\cr
getRootPath

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
as.character, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, print


 }

 \section{Requirements}{
   This class requires that an aroma probe sequence file and aroma
   match scores file is available for the chip type.
 }

 \author{
   Mark Robinson, WEHI.
 }

 \references{
   [1] Johnson WE, Li W, Meyer CA, Gottardo R, Carroll JS, Brown M, Liu XS.
     \emph{Model-based analysis of tiling-arrays for ChIP-chip}, PNAS, 2006.
 }
\keyword{classes}
