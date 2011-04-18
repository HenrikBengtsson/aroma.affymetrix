%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ScaleNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ScaleNormalization}
\docType{class}
\alias{ScaleNormalization}


 \title{The ScaleNormalization class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ScaleNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{ScaleNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ScaleNormalization}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr



  This class represents a normalization function that transforms the
  probe-level signals towards the same scale.
 }

 \usage{ScaleNormalization(dataSet=NULL, ..., targetAvg=4400, subsetToUpdate=NULL, typesToUpdate=NULL, subsetToAvg="-XY", typesToAvg=typesToUpdate, shift=0)}

 \arguments{
   \item{dataSet}{\code{\link{AffymetrixCelSet}} to be normalized.}
   \item{...}{Arguments passed to the constructor of
     \code{\link{ProbeLevelTransform}}.}
   \item{targetAvg}{A \code{\link[base]{numeric}} value.}
   \item{subsetToUpdate}{The probes to be updated.
     If \code{\link[base]{NULL}}, all probes are updated.}
   \item{typesToUpdate}{Types of probes to be updated.}
   \item{subsetToAvg}{The probes to calculate average signal over.
     If a single \code{\link[base]{numeric}} in (0,1), then this
     fraction of all probes will be used.
     If \code{\link[base]{NULL}}, all probes are considered.}
   \item{typesToAvg}{Types of probes to be used when calculating the
     average signal.
     If \code{"pm"} and \code{"mm"} only perfect-match and mismatch
     probes are used, respectively. If \code{"pmmm"} both types are used.
   }
   \item{shift}{Optional amount of shift if data before fitting/normalizing.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.ScaleNormalization]{process}} \tab Normalizes the data set.\cr
}


 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputDataSetOLD20090509, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
