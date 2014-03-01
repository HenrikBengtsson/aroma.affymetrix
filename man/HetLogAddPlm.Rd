%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  HetLogAddPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{HetLogAddPlm}
\docType{class}
\alias{HetLogAddPlm}


\title{The HetLogAddPlm class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class HetLogAddPlm}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{HetLogAddPlm}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{HetLogAddCnPlm}}, \emph{\link[aroma.affymetrix]{HetLogAddSnpPlm}}\cr

public abstract static class \bold{HetLogAddPlm}\cr
extends \emph{\link[aroma.affymetrix]{RmaPlm}}\cr



 This class represents a log-additive model similar to the one described
 in Irizarry et al (2003), except that the errors may have different
 variances for different probes.
}

\usage{
HetLogAddPlm(...)
}

\arguments{
  \item{...}{Arguments passed to \code{\link{RmaPlm}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


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

\author{Henrik Bengtsson}

\seealso{
 \code{\link{RmaPlm}}.
}
\keyword{classes}
