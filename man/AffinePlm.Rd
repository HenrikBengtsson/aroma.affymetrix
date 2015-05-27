%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffinePlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AffinePlm}
\docType{class}
\alias{AffinePlm}

\encoding{latin1}

\title{The AffinePlm class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class AffinePlm}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AffinePlm}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffineCnPlm}}, \emph{\link[aroma.affymetrix]{AffineSnpPlm}}\cr

public abstract static class \bold{AffinePlm}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelModel}}\cr



 This class represents affine model in Bengtsson \& Hossjer (2006).
}

\usage{
AffinePlm(..., background=TRUE)
}

\arguments{
  \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
  \item{background}{If \code{\link[base:logical]{TRUE}}, background is estimate for each unit group,
    otherwise not. That is, if \code{\link[base:logical]{FALSE}}, a \emph{linear} (proportional)
    model without offset is fitted, resulting in very similar results as
    obtained by the \code{\link{MbeiPlm}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getProbeAffinityFile} \tab  -\cr
}


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
  For a single unit group, the affine model is:

   \deqn{y_{ik} = a + \theta_i \phi_k + \varepsilon_{ik}}

  where \eqn{a} is an offset common to all probe signals,
  \eqn{\theta_i} are the chip effects for arrays \eqn{i=1,...,I},
  and \eqn{\phi_k} are the probe affinities for probes \eqn{k=1,...,K}.
  The \eqn{\varepsilon_{ik}} are zero-mean noise with equal variance.
  The model is constrained such that \eqn{\prod_k \phi_k = 1}.

  Note that with the additional constraint \eqn{a=0} (see arguments above),
  the above model is very similar to \code{\link{MbeiPlm}}.  The differences in
  parameter estimates is due to difference is assumptions about the
  error structure, which in turn affects how the model is estimated.
}

\author{Henrik Bengtsson}

\references{
  Bengtsson \& Hossjer (2006). \cr
}
\keyword{classes}
