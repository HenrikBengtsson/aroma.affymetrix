%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GcContentNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GcContentNormalization}
\docType{class}
\alias{GcContentNormalization}


 \title{The GcContentNormalization class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GcContentNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectTransform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{GcContentNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{GcContentNormalization}\cr
extends \emph{\link[aroma.affymetrix]{ChipEffectTransform}}\cr


 }

 \usage{GcContentNormalization(dataSet=NULL, ..., targetFunction=NULL, subsetToFit=NULL)}

 \arguments{
   \item{dataSet}{A \code{\link{CnChipEffectSet}}.}
   \item{...}{Additional arguments passed to the constructor of
     \code{\link{ChipEffectTransform}}.}
   \item{targetFunction}{A \code{\link[base]{function}}. The target function to which all arrays
     should be normalized to.}
   \item{subsetToFit}{The units from which the normalization curve should
     be estimated.  If \code{\link[base]{NULL}}, all are considered.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getOutputDataSet00} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:process.GcContentNormalization]{process}} \tab Normalizes the data set.\cr
}


 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputDataSetOLD20090509, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \section{Requirements}{
   This class requires an Aroma unit GC-content (UGC) file.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
