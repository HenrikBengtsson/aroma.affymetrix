%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ExonRmaPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ExonRmaPlm}
\docType{class}
\alias{ExonRmaPlm}


 \title{The ExonRmaPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ExonRmaPlm}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{ExonRmaPlm}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ExonRmaPlm}\cr
extends \link[aroma.affymetrix]{RmaPlm}\cr



  This class represents the log-additive model part of the Robust Multichip
  Analysis (RMA) method described in Irizarry et al (2003), as implemented
  for exon arrays.  The model may be fitted with exons merged into
  transcripts (all probes fitted together) or on an individual exon basis
  (probes within an exon treated as a group, but exons fitted separately).
 }

 \usage{ExonRmaPlm(..., mergeGroups=TRUE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{RmaPlm}}.}
   \item{mergeGroups}{A \code{\link[base]{logical}} flag specifying whether to merge exons
      into transcripts.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{calculateWeights} \tab  -\cr
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{getChipEffects} \tab  -\cr
 \tab \code{getChipEffectSet} \tab  -\cr
 \tab \code{getProbeAffinities} \tab  -\cr
 \tab \code{getProbeAffinityFile} \tab  -\cr
 \tab \code{setMergeGroups} \tab  -\cr
}


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
    \code{\link{RmaPlm}}.
 }

 \author{Ken Simpson (ksimpson[at]wehi.edu.au).}

 \references{
  Irizarry et al. \emph{Summaries of Affymetrix GeneChip probe level data}.
  NAR, 2003, 31, e15.\cr
 }
\keyword{classes}
