%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpChipEffectGroupMerge.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SnpChipEffectGroupMerge}
\docType{class}
\alias{SnpChipEffectGroupMerge}


 \title{The SnpChipEffectGroupMerge class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class SnpChipEffectGroupMerge}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectTransform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectGroupMerge}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{SnpChipEffectGroupMerge}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{SnpChipEffectGroupMerge}\cr
extends \link[aroma.affymetrix]{ChipEffectGroupMerge}\cr



  This class represents a method that merges SNP chip effects across groups
  unit by unit.
 }

 \usage{SnpChipEffectGroupMerge(..., mergeStrands=FALSE, mean=c("arithmetic", "geometric"))}

 \arguments{
   \item{...}{Arguments passed to the constructor of
        \code{\link{ChipEffectGroupMerge}}.}
   \item{mergeStrands}{If \code{\link[base:logical]{TRUE}}, group strands are merged.}
   \item{mean}{A \code{\link[base]{character}} string specifying what type of averaging
        should be applied.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getAsteriskTags} \tab  -\cr
 \tab \code{getMergeFunction} \tab  -\cr
}


 \bold{Methods inherited from ChipEffectGroupMerge}:\cr
getCdf, getMergeFunction, process, setMergeFunction

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputDataSetOLD20090509, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{classes}
\keyword{internal}
