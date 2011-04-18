%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipQuantileNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{DChipQuantileNormalization}
\docType{class}
\alias{DChipQuantileNormalization}


 \title{The DChipQuantileNormalization class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class DChipQuantileNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{QuantileNormalization}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{DChipQuantileNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{DChipQuantileNormalization}\cr
extends \link[aroma.affymetrix]{QuantileNormalization}\cr



  This class represents a special \code{\link{QuantileNormalization}}
  using smooth-splines.
 }

 \usage{DChipQuantileNormalization(..., robust=FALSE)}

 \arguments{
   \item{...}{Arguments passed to the constructor of
       \code{\link{QuantileNormalization}}.}
   \item{robust}{If \code{\link[base:logical]{TRUE}}, the normalization function is estimated
       robustly, otherwise not.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.DChipQuantileNormalization]{process}} \tab Normalizes the data set.\cr
}


 \bold{Methods inherited from QuantileNormalization}:\cr
clearCache, findTargetDistributionFile, process

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputDataSetOLD20090509, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \details{
   This normalization method implements the two-pass algorithm described
   in Bengtsson et al. (2008).
 }

 \references{
   [1] H. Bengtsson, R. Irizarry, B. Carvalho, & T.P. Speed.
       Estimation and assessment of raw copy numbers at the single
       locus level, Bioinformatics, 2008.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
