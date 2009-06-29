%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  UnitTypeScaleNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{UnitTypeScaleNormalization}
\docType{class}
\alias{UnitTypeScaleNormalization}


 \title{The UnitTypeScaleNormalization class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class UnitTypeScaleNormalization}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{UnitTypeScaleNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{UnitTypeScaleNormalization}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelTransform3}}\cr



  This class represents a normalization function that transforms the
  probe signals such that each unit type gets the same average.
 }

 \usage{UnitTypeScaleNormalization(..., targetAvg=4400)}

 \arguments{
   \item{...}{Arguments passed to the constructor of
     \code{\link{ProbeLevelTransform3}}.}
   \item{targetAvg}{A \code{\link[base]{numeric}} value.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.UnitTypeScaleNormalization]{process}} \tab Normalizes the data set.\cr
}


 \bold{Methods inherited from ProbeLevelTransform3}:\cr
clearCache, getCellsTo, getCellsToFit, getCellsToUpdate, getUnitsTo, getUnitsToFit, getUnitsToUpdate, writeSignals

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputDataSetOLD20090509, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \examples{\dontrun{
 }}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
