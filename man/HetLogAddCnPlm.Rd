%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  HetLogAddCnPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{HetLogAddCnPlm}
\docType{class}
\alias{HetLogAddCnPlm}


\title{The HetLogAddCnPlm class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class HetLogAddCnPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{MultiArrayUnitModel}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{RmaPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{HetLogAddPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{HetLogAddSnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{CnPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{HetLogAddCnPlm}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{HetLogAddCnPlm}\cr
extends \link[aroma.affymetrix]{CnPlm}\cr


}

\usage{
HetLogAddCnPlm(..., combineAlleles=FALSE)
}

\arguments{
  \item{...}{Arguments passed to \code{\link{HetLogAddCnPlm}}.}
  \item{combineAlleles}{If \code{\link[base:logical]{FALSE}}, allele A and allele B are treated
     seperately, otherwise together.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from CnPlm}:\cr
getCellIndices, getChipEffectSet, getCombineAlleles, getParameters, getProbeAffinityFile, setCombineAlleles

 \bold{Methods inherited from SnpPlm}:\cr
getCellIndices, getChipEffectSet, getMergeStrands, getParameters, getProbeAffinityFile, setMergeStrands

 \bold{Methods inherited from HetLogAddSnpPlm}:\cr
getAsteriskTags

 \bold{Methods inherited from HetLogAddPlm}:\cr
getAsteriskTags, getFitUnitGroupFunction

 \bold{Methods inherited from RmaPlm}:\cr
getAsteriskTags, getCalculateResidualsFunction, getParameters, getRlmFitFunctions

 \bold{Methods inherited from ProbeLevelModel}:\cr
calculateResidualSet, calculateWeights, fit, getAsteriskTags, getCalculateResidualsFunction, getChipEffectSet, getProbeAffinityFile, getResidualSet, getRootPath, getWeightsSet

 \bold{Methods inherited from MultiArrayUnitModel}:\cr
getListOfPriors, setListOfPriors, validate

 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getAsteriskTags, getFitSingleCellUnitFunction, getParameters

 \bold{Methods inherited from Model}:\cr
as.character, fit, getAlias, getAsteriskTags, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, asThis


}

\section{Model}{
  TO DO.
}

\author{Henrik Bengtsson}
\keyword{classes}
