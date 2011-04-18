%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RmaCnPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{RmaCnPlm}
\docType{class}
\alias{RmaCnPlm}


 \title{The RmaCnPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class RmaCnPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{MultiArrayUnitModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{RmaPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{RmaSnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{CnPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{RmaCnPlm}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{RmaCnPlm}\cr
extends \link[aroma.affymetrix]{CnPlm}\cr


 }

 \usage{RmaCnPlm(..., combineAlleles=FALSE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{RmaSnpPlm}}.}
   \item{combineAlleles}{If \code{\link[base:logical]{FALSE}}, allele A and allele B are treated
      seperately, otherwise together.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from CnPlm}:\cr
getCellIndices, getChipEffectSet, getCombineAlleles, getProbeAffinityFile, setCombineAlleles

 \bold{Methods inherited from SnpPlm}:\cr
getCellIndices, getChipEffectSet, getMergeStrands, getProbeAffinityFile, setMergeStrands

 \bold{Methods inherited from RmaSnpPlm}:\cr
getAsteriskTags

 \bold{Methods inherited from RmaPlm}:\cr
getAsteriskTags, getCalculateResidualsFunction, getRlmFitFunctions

 \bold{Methods inherited from ProbeLevelModel}:\cr
calculateResidualSet, calculateWeights, fit, getAsteriskTags, getCalculateResidualsFunction, getChipEffectSet, getProbeAffinityFile, getResidualSet, getWeightsSet

 \bold{Methods inherited from MultiArrayUnitModel}:\cr
getListOfPriors, setListOfPriors, validate

 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getAsteriskTags, getFitSingleCellUnitFunction

 \bold{Methods inherited from Model}:\cr
fit, getAlias, getAsteriskTags, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \section{Model}{
   TO DO.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
