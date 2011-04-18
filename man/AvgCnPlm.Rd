%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AvgCnPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AvgCnPlm}
\docType{class}
\alias{AvgCnPlm}


 \title{The AvgCnPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AvgCnPlm}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AvgPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AvgSnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{CnPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{AvgCnPlm}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AvgCnPlm}\cr
extends \link[aroma.affymetrix]{CnPlm}\cr


 }

 \usage{AvgCnPlm(..., combineAlleles=FALSE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AvgSnpPlm}}.}
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

 \bold{Methods inherited from AvgSnpPlm}:\cr
getAsteriskTags

 \bold{Methods inherited from AvgPlm}:\cr
getAsteriskTags, getCalculateResidualsFunction, validate

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

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
