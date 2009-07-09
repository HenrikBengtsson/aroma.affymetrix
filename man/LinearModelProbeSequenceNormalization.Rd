%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  LinearModelProbeSequenceNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{LinearModelProbeSequenceNormalization}
\docType{class}
\alias{LinearModelProbeSequenceNormalization}


 \title{The LinearModelProbeSequenceNormalization class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class LinearModelProbeSequenceNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform3}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AbstractProbeSequenceNormalization}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{LinearModelProbeSequenceNormalization}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{BasePositionNormalization}\cr

public abstract static class \bold{LinearModelProbeSequenceNormalization}\cr
extends \emph{\link[aroma.affymetrix]{AbstractProbeSequenceNormalization}}\cr



  This abstract class represents a normalization method that corrects
  for systematic effects in the probe intensities due to probe-sequence
  dependent effects that can be modelled using a linear model.
 }

 \usage{LinearModelProbeSequenceNormalization(...)}

 \arguments{
   \item{...}{Arguments passed to the constructor of
     \code{\link{AbstractProbeSequenceNormalization}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from AbstractProbeSequenceNormalization}:\cr
fitOne, getAromaCellSequenceFile, getTargetFile, indexOfMissingSequences, predictOne, process

 \bold{Methods inherited from ProbeLevelTransform3}:\cr
clearCache, getCellsTo, getCellsToFit, getCellsToUpdate, getUnitsTo, getUnitsToFit, getUnitsToUpdate, writeSignals

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputDataSetOLD20090509, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \section{Requirements}{
   This class requires that an aroma probe sequence file is available
   for the chip type.
 }

 \section{Memory usage}{
  The model fitting methods of this class are bounded in memory.
  This is done by first building up the normal equations incrementally
  in chunks of cells.  The generation of normal equations is otherwise
  the step that consumes the most memory.
  When the normal equations are available, the \code{\link[base]{solve}}()
  method is used to solve the equations.  Note that this algorithm is
  still exact.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
